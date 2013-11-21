// __BEGIN_LICENSE__
//  Copyright (c) 2009-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NGT platform is licensed under the Apache License, Version 2.0 (the
//  "License"); you may not use this file except in compliance with the
//  License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


/// \file lronacSolverModel.cc
///

#include <iostream>

#include <iTime.h> // Isis time class

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp> // for null_deleter

#include <vw/InterestPoint.h>
#include <vw/Image/MaskViews.h>
//#include <boost/accumulators/accumulators.hpp>
//#include <boost/accumulators/statistics.hpp>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Stereo/PreFilter.h>
//#include <vw/Stereo/CorrelationView.h>
//#include <vw/Stereo/CostFunctions.h>
//#include <vw/Stereo/DisparityMap.h>

#include <vw/Stereo/Correlate.h>

//#include <asp/Core/DemDisparity.h>
//#include <asp/Core/LocalHomography.h>

#include <vw/Stereo/StereoModel.h>
#include <asp/IsisIO/IsisCameraModel.h> //::point_to_pixel>
#include <vw/Math/LevenbergMarquardt.h>
#include <vw/Math/EulerAngles.h>

#include <asp/Core/IntegralAutoGainDetector.h>
#include <asp/IsisIO/IsisInterfaceLineScan.h>

#include <stereo.h>

#include <lronacSolverSupport.h>

using namespace vw;
using namespace vw::stereo;
using namespace asp::isis;
using namespace asp;
using std::endl;
using std::setprecision;
using std::setw;


//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// Ceres model code




#include "ceres/ceres.h"
#include "glog/logging.h"


/// Class for solving for the rotation between two LRONAC cameras
class LrocPairModel : public vw::math::LeastSquaresModelBase<LrocPairModel>
{
public:  // Definitions ------------------------------------------------------------------------------

/// * Defines a result_type that is the type returned by
///   evaluating the functor.  Typically Vector<float> or
///   Vector<double>
typedef Vector<double> result_type;

/// * Defines a domain_type that is the type of the search
///   space.  Often a Vector<foo>, but can reflect other
///   topologies if needed.
typedef Vector<double> domain_type;


/// * Defines a jacobian_type corresponding to the space of
///   jacobian matrices.  Typically Matrix<foo>.
typedef Matrix<double> jacobian_type;

private: // Variables ---------------------------------------------------------------------------

// Camera models
asp::isis::IsisInterfaceLineScanRot _leftCameraModel;
asp::isis::IsisInterfaceLineScanRot _rightCameraModel;

bool _solveWorldFrame;
bool _includePosition; // Currently only works if solveWorldFrame is true!

mutable vw::camera::AdjustedCameraModel _rightCameraRotatedModel;

//TODO: Modify this class so that the rays from this task can be used!
vw::stereo::StereoModel         _stereoModel;

// Observation records
Vector<double> _leftRows;
Vector<double> _rightRows;
Vector<double> _leftCols;
Vector<double> _rightCols;

public: // Functions  -----------------------------------------------------------------------------------

/// Constructor performs initialization
LrocPairModel(const std::string &leftCubePath,  const std::string &rightCubePath, 
              const bool solveWorldFrame=false, const bool includePosition=false)
  : _leftCameraModel(leftCubePath), _rightCameraModel(rightCubePath), 
    _solveWorldFrame(solveWorldFrame), _includePosition(includePosition),
    _rightCameraRotatedModel(boost::shared_ptr<vw::camera::CameraModel>(&_rightCameraModel, boost::serialization::null_deleter())),
    _stereoModel(&_leftCameraModel, &_rightCameraModel) // Only the right camera is rotated
{
  // Both camera models are loaded from file on initialization
  printf("Done constructing LROC model\n");
}

/// Test functions to get a single pixel vector
Vector3 getLeftVector (const Vector2& pixel) {return vw::math::normalize(_leftCameraModel.pixel_to_vector (pixel));}
/*  Vector3 getRightVector(const Vector2& pixel) 
{
  if (_solveWorldFrame)
    return vw::math::normalize(_rightCameraRotatedModel.pixel_to_vector(pixel));
  else // Camera frame
    return vw::math::normalize(_rightCameraModel.pixel_to_vector(pixel));
}
*/

Vector2 getRightPixelRot(const Vector3& point, const Vector3& rot, const Vector3& offset=Vector3())
{
  if (_solveWorldFrame)
  {
    _rightCameraRotatedModel.set_axis_angle_rotation(rot);
    _rightCameraRotatedModel.set_translation(offset);
    return vw::math::normalize(_rightCameraRotatedModel.point_to_pixel(point));
  }
  else // Camera frame
    return vw::math::normalize(_rightCameraModel.point_to_pixel_rotated(point, rot));
}


/// Given the pixel pair observations, compute the initial state estimate.
/// - This also loads the observation vectors and returns a packed version of them.
bool getInitialStateEstimate(const Vector<double> &leftRows,  const Vector<double> &leftCols,
                              const Vector<double> &rightRows, const Vector<double> &rightCols,
                              Vector<double> &stateEstimate, Vector<double> &packedObsVector,
                              const std::vector<double> &inputState = std::vector<double>())
{ 
  //TODO: Verify vector sizes match
  
  // Record observation vectors to class variables
  _leftRows  = leftRows;
  _rightRows = rightRows;
  _leftCols  = leftCols;
  _rightCols = rightCols;    // Maybe don't need these

  /* Input parameter set:
    rotationOffsetX (additional rotations applied to RE camera)
    rotationOffsetY
    rotationOffsetZ
    [positionOffsetX // Optional with certain parameters
      positionOffsetY
      positionOffsetZ
    ]
    x               [Repeated for every correspondence point]
    y
    z
  */
  
  // Determine the number of state elements
  const size_t numPoints        = leftRows.size();
        size_t numFixedElements = 3;
  if (_includePosition) 
    numFixedElements = 6;
  const size_t numStateElements = numPoints*3 + numFixedElements;
  stateEstimate.set_size(numStateElements);
  packedObsVector.set_size(numPoints*4);

  // Rotation values start at zero
  stateEstimate[0] = 0;
  stateEstimate[1] = 0;
  stateEstimate[2] = 0;

  if (_includePosition) // Also init the position to zero
  {
    stateEstimate[3] = 0;
    stateEstimate[4] = 0;
    stateEstimate[5] = 0;
  }

  // Unless an input state was given start at zero rotation
  Vector3 zeroVec(0, 0, 0); 
  Vector3 rotVec(stateEstimate[0], stateEstimate[1], stateEstimate[2]);

  // If an initial state was provided, use those values instead of zeroes.
  if (inputState.size() > 0)
  {
    printf("Calculating initial points using input state...\n");
    for (size_t i=0; i<inputState.size(); ++i)
      stateEstimate[i] = inputState[i];

    if (_solveWorldFrame) // Update _rightCameraRotatedModel
    {
      rotVec = Vector3(stateEstimate[0], stateEstimate[1], stateEstimate[2]);
      _rightCameraRotatedModel.set_axis_angle_rotation(rotVec);
      if (_includePosition) // Also update translation
      {
        Vector3 offsetVec(stateEstimate[3], stateEstimate[4], stateEstimate[5]);
        _rightCameraRotatedModel.set_translation(offsetVec);
      }
    }
    else //TODO: Support this case one day
    {
      printf("Initial state with local transform is not currently supported!\n");
      return false;
    }
  } 
  else // Default starting state
  {
    _rightCameraRotatedModel.set_axis_angle_rotation(zeroVec);
    _rightCameraRotatedModel.set_translation(zeroVec);
  }

  //DEBUG
  // Set up georeference class with default moon datum
  vw::cartography::Datum datum("D_MOON");

  vw::cartography::GeoReference leftGeoRef(datum);
  vw::cartography::GeoReference rightGeoRef(datum);

  

  printf("Setting up state estimate\n");

  // For each input point pair
  Vector3 pointLoc, lastPointLoc;
  for (size_t i=0; i<numPoints; ++i)
  {
    // Set up point pair
    Vector2 leftPixel (leftCols [i], leftRows [i]);
    Vector2 rightPixel(rightCols[i], rightRows[i]);
    
    // Compute the intersection location
    double triangulationError;
    //pointLoc = _stereoModel(leftPixel, rightPixel, triangulationError); // Not working!
    
    Vector3 leftCamCenter  = _leftCameraModel.camera_center(leftPixel);
    Vector3 rightCamCenter = _rightCameraRotatedModel.camera_center(rightPixel);
    
    Vector3            leftVec      = vw::math::normalize(_leftCameraModel.pixel_to_vector (leftPixel));
    Vector3            rightVec     = vw::math::normalize(_rightCameraRotatedModel.pixel_to_vector(rightPixel));
    double             vectorAngle  = acos(vw::math::dot_prod(leftVec, rightVec));
    
    // Quaternion code did not seem to be working
    //Quaternion<double> leftCamPose  = vw::math::normalize(_leftCameraModel.camera_pose (leftPixel));
    //Quaternion<double> rightCamPose = vw::math::normalize(_rightCameraModel.camera_pose(rightPixel));
    //Quaternion<double> conjLeftPose = vw::math::conj(leftCamPose);
    //Quaternion<double> invLeftPose  = vw::math::inverse(leftCamPose);
    //Quaternion<double> rotDiff      = rightCamPose*invLeftPose;
    //Vector3 axis; double angle;
    //rotDiff.axis_angle(axis, angle);

    //std::cout << "leftCamVector "  << leftVec      << std::endl;
    //std::cout << "rightCamVector " << rightVec     << std::endl;      
    //std::cout << "camVectorDiff "  << leftVec-rightVec << std::endl;  
    //std::cout << "vectorAngle (degrees) "    << vectorAngle*180/3.15159  << std::endl;
    //std::cout << "leftCamPose "    << leftCamPose  << std::endl;
    //std::cout << "rightCamPose "   << rightCamPose << std::endl;
    //std::cout << "conjLeftPose "   << conjLeftPose << std::endl;
    //std::cout << "invLeftPose "    << invLeftPose  << std::endl;
    //std::cout << "rotDiff "        << rotDiff      << std::endl;
    
Vector3 v12 = cross_prod(leftVec, rightVec);
Vector3 v1  = cross_prod(v12,     leftVec);
Vector3 v2  = cross_prod(v12,     rightVec);

Vector3 closestPoint1 = leftCamCenter  + dot_prod(v2, rightCamCenter-leftCamCenter )/dot_prod(v2, leftVec )*leftVec;
Vector3 closestPoint2 = rightCamCenter + dot_prod(v1, leftCamCenter -rightCamCenter)/dot_prod(v1, rightVec)*rightVec;

const double MOON_RADIUS_M = 1737400;
if (inputState.empty()) // Override the stereo results with ground intersections!
{
  // Things don't work if this is done with an initial state
  if (!raySphereIntersect(leftCamCenter,  leftVec,  MOON_RADIUS_M, closestPoint1))
    printf("Failed to intersect moon with left!\n");
  if (!raySphereIntersect(rightCamCenter, rightVec, MOON_RADIUS_M, closestPoint2))
    printf("Failed to intersect moon with right!\n");
}
else
{
  if ((i % 1000) == 0) // Test output
    printf("%d\n", i);
}

Vector3 midPoint = 0.5 * (closestPoint1 + closestPoint2);
pointLoc = midPoint; // HIJACK TRIANGULATION CALCULATIONS!
    

    
    
Vector3 errorVec = closestPoint1 - closestPoint2;
triangulationError = vw::math::norm_2(errorVec);
        
    //printf("\n");
    //std::cout << "closestPoint1   = " << closestPoint1 <<", radius = " << vw::math::norm_2(closestPoint1)/1000.0 << std::endl;
    //std::cout << "closestPoint2   = " << closestPoint2 <<", radius = " << vw::math::norm_2(closestPoint2)/1000.0 << std::endl;
    //std::cout << "errorVec        = " << errorVec      << std::endl;
    //std::cout << "projection dist (km)= " << vw::math::norm_2(closestPoint1 - leftCamCenter)/1000.0 << std::endl;
        
    
    //std::cout.precision(16);
    //std::cout << "Left  camera center = " << leftCamCenter [0] <<", "<< leftCamCenter [1] <<", "<< leftCamCenter [2]<<", radius = " << vw::math::norm_2(leftCamCenter)/1000.0 << std::endl;
    //std::cout << "Right camera center = " << rightCamCenter[0] <<", "<< rightCamCenter[1] <<", "<< rightCamCenter[2]<<", radius = " << vw::math::norm_2(rightCamCenter)/1000.0 << std::endl;
    //std::cout << "Center diff         = " << leftCamCenter[0]-rightCamCenter[0] <<", "<< leftCamCenter[1]-rightCamCenter[1] <<", "<< leftCamCenter[2]-rightCamCenter[2] << std::endl;
    //printf("Angle between cameras = %lf\n", vectorAngle);
    //printf("Center abs diff = %lf\n", sqrt( pow(leftCamCenter[0]-rightCamCenter[0], 2) + pow(leftCamCenter[1]-rightCamCenter[1], 2) + pow(leftCamCenter[2]-rightCamCenter[2], 2) ));
    
    
    //std::cout << "Left pixel  = " << leftPixel << " Right pixel = " << rightPixel << std::endl;
    //std::cout << "Initial point " << i << " = " << pointLoc << " Triangulation error = " << triangulationError << std::endl;


    // DEBUG code to project the computed intersection back into the cameras
    Vector2 rightProjection, leftProjection;
    if (_solveWorldFrame)
      rightProjection = _rightCameraRotatedModel.point_to_pixel(pointLoc);
    else // Camera frame
      rightProjection = _rightCameraModel.point_to_pixel_rotated(pointLoc, rotVec);
    leftProjection  = _leftCameraModel.point_to_pixel(pointLoc);

    //std::cout << "**Left: "  << leftProjection  << " - " << leftPixel  << std::endl;
    //std::cout << "**Right: " << rightProjection << " - " << rightPixel << std::endl;
    
    // Sanity check
    const double intersectionRadius = vw::math::norm_2(pointLoc);
    const double camRadius          = vw::math::norm_2(rightCamCenter);
    //std::cout << "intersection radius = " << intersectionRadius << std::endl;
    
    // Convert from GCC to GDC
    //Vector3 gdcCoord = datum.cartesian_to_geodetic(midPoint);
    //std::cout << "elevation above ellipsoid = " << gdcCoord[2] << std::endl;
    
    // Try dropping the elevation down to the datum level
    //gdcCoord[2] = 0.0;
    //pointLoc = datum.geodetic_to_cartesian(gdcCoord);
    
    if (camRadius < intersectionRadius)
    {
      printf("Warning: Point %lu, reverse intersection!  Using previous point location as an estimate.\n", i);
      printf("%lf < %lf\n", camRadius, intersectionRadius);
      if (i == 0) // No previous point to copy
        pointLoc = leftCamCenter + (leftVec*37000); // Extend 37km from left camera
      else
        pointLoc = lastPointLoc; // For now our best guess to a failed projection is the previous point projection (at least it is in front of the camera!)
    } 
    else // Point estimated successfully
      lastPointLoc = pointLoc;
    
    // Record the x/y/z value for this point
    stateEstimate[numFixedElements + i*3 + 0] = pointLoc[0];
    stateEstimate[numFixedElements + i*3 + 1] = pointLoc[1];
    stateEstimate[numFixedElements + i*3 + 2] = pointLoc[2];
    
    // Build the packed observation vector
    packedObsVector[i*4 + 0] = leftCols [i];
    packedObsVector[i*4 + 1] = leftRows [i];
    packedObsVector[i*4 + 2] = rightCols[i];
    packedObsVector[i*4 + 3] = rightRows[i];
    
  } // End loop through points
  
  return true;
} // end getInitialStateEstimate()


/// Generates a vector containing the stereo intersection distance for each pixel pair
/// - This gives a good indication of how accurate the current correction is
std::vector<double> computeError(domain_type const& x)
{    
  // Determine the number of state elements
  const size_t numPoints = _leftRows.size();
  std::vector<double> errorVector(numPoints);

  Vector3 rotVec   (x[0], x[1], x[2]);
  Vector3 offsetVec(x[3], x[4], x[5]); // Only used if _includePosition = true

  int firstPointState = 3;
  if (_includePosition) 
  {
    firstPointState = 6;
    _rightCameraRotatedModel.set_translation(offsetVec);
  }
  if (_solveWorldFrame)
    _rightCameraRotatedModel.set_axis_angle_rotation(rotVec);

      // For each input point pair
  Vector3 pointLoc, lastPointLoc;
  for (size_t i=0; i<numPoints; ++i)
  {
    // Create a point object
    Vector3 thisPoint(x[firstPointState + i*3 + 0],
                      x[firstPointState + i*3 + 1],
                      x[firstPointState + i*3 + 2]);      
    
    //std::cout << "thisPoint = " << thisPoint << std::endl;
                      
    // Project to pixel locations                  
    Vector2 leftProjection  = _leftCameraModel.point_to_pixel(thisPoint);
    Vector2 rightProjection;
    if (_solveWorldFrame)
      rightProjection = _rightCameraRotatedModel.point_to_pixel(thisPoint); // This treats the angles as a quaternion vector!
    else // Camera frame
      rightProjection = _rightCameraModel.point_to_pixel_rotated(thisPoint, rotVec); // This treats the angles as euler angles!
    
                      
    // Compare to observed pixel locations
    Vector2 leftObsPoint (_leftCols [i], _leftRows [i]);
    Vector2 rightObsPoint(_rightCols[i], _rightRows[i]);
    
    Vector2 leftDiff  = leftProjection  - leftObsPoint;
    Vector2 rightDiff = rightProjection - rightObsPoint;
    //std::cout << "Left: "  << leftProjection  << " - " << leftObsPoint  << std::endl;
    //std::cout << "Right: " << rightProjection << " - " << rightObsPoint << std::endl;
    
    double leftError  = vw::math::norm_2(leftDiff );
    double rightError = vw::math::norm_2(rightDiff);
    
    // Returned value is average of left and right position error
    errorVector[i] = (leftError + rightError) / 2.0;
  } // End loop through points
  
  return errorVector;
    
} // end computeError()

/// For a given point, compute the left camera observation.
/// - Returns false if the point is not visible.
bool getLeftObservation(const double* const pointParams, double *observation) const
{
  // Create a point object
  Vector3 thisPoint(pointParams[0], pointParams[1], pointParams[2]);
  Vector2 leftProjection; 
  try // Project the point into the camera
  {
    leftProjection  = _leftCameraModel.point_to_pixel(thisPoint);
  }
  catch(...) // Handle errors
  {
    std::cout << "Warning: Failed to project location: " << thisPoint << ", using previous intersection location." << std::endl;
    return false;
  }
  
  // Load the projected pixels into the output observation vector
  observation[0] = leftProjection[0]; // x
  observation[1] = leftProjection[1]; // y

  return true;
}

/// For a given point, compute the right camera observation.
/// - Returns false if the point is not visible.
bool getRightObservation(const double* const cameraParams, const double* const pointParams, double *observation)
{
  // This function returns an error vector for a given set of parameters

  // Apply the rotations from the state vector to the right LROC camera model
  Vector3 rotVec(cameraParams[0], cameraParams[1], cameraParams[2]);

  if (_includePosition) 
  {
    Vector3 offsetVec(cameraParams[3], cameraParams[4], cameraParams[5]); 
    _rightCameraRotatedModel.set_translation(offsetVec);
  }
  if (_solveWorldFrame)
    _rightCameraRotatedModel.set_axis_angle_rotation(rotVec);

  // Create a point object
  Vector3 thisPoint(pointParams[0], pointParams[1], pointParams[2]);
  Vector2 rightProjection;
  try // Project the point into the camera
  {
    if (_solveWorldFrame)
      rightProjection = _rightCameraRotatedModel.point_to_pixel(thisPoint);
    else // Camera frame
      rightProjection = _rightCameraModel.point_to_pixel_rotated(thisPoint, rotVec);
  }
  catch(...)
  {
    std::cout << "Warning: Failed to project location: " << thisPoint << ", using previous intersection location." << std::endl;
    return false;
  }

  // Load the projected pixels into the output obserVation vector
  observation[0] = rightProjection[0];
  observation[1] = rightProjection[1];

  return true;
}


/// * Defines a method: result_type operator()( domain_type const& x ) const;
///   that evaluates the model function at the given point.
result_type operator()( domain_type const& x ) const
{
  // This function returns an error vector for a given set of parameters

  // Apply the rotations from the state vector to the right LROC camera model
  Vector3 rotVec   (x[0], x[1], x[2]);
  Vector3 offsetVec(x[3], x[4], x[5]); // Only used if _includePosition = true

  int firstPointState = 3;
  if (_includePosition) 
  {
    firstPointState = 6;
    _rightCameraRotatedModel.set_translation(offsetVec);
  }
  if (_solveWorldFrame)
    _rightCameraRotatedModel.set_axis_angle_rotation(rotVec);
  
  const double rad2deg = 180.0 / M_PI;
  //printf("Trying rotation %lf, %lf, %lf\n", x[0]*rad2deg, x[1]*rad2deg, x[2]*rad2deg);
  
  //std::cout << "Left  camera center = " << _leftCameraModel.camera_center()  << std::endl;
  //std::cout << "Right camera center = " << _rightCameraModel.camera_center() << std::endl;
  
  // Set up output vector
  size_t numPoints = _leftRows.size();
  Vector<double> obsVec(numPoints*4);
  
  // Compute expected obseration value at each pixel
  for (size_t i=0; i<numPoints; ++i)
  {
    // Create a point object
    Vector3 thisPoint(x[firstPointState + i*3 + 0],
                      x[firstPointState + i*3 + 1],
                      x[firstPointState + i*3 + 2]);
//      std::cout << "This point = " << thisPoint << std::endl;

    // Project the point into both cameras
    Vector2 leftProjection, rightProjection;
    try
    {
      leftProjection  = _leftCameraModel.point_to_pixel(thisPoint);
      if (_solveWorldFrame)
        rightProjection = _rightCameraRotatedModel.point_to_pixel(thisPoint);
      else // Camera frame
        rightProjection = _rightCameraModel.point_to_pixel_rotated(thisPoint, rotVec);
    }
    catch(...)
    {
      std::cout << "Warning: Failed to project location " << i << ": " << thisPoint << ", using previous intersection location." << std::endl;
      
      //TODO: Better solution than just using the last location!
      if (i == 0)
        return false;
      obsVec[4*i + 0] = obsVec[4*(i-1) + 0];
      obsVec[4*i + 1] = obsVec[4*(i-1) + 1];
      obsVec[4*i + 2] = obsVec[4*(i-1) + 2];
      obsVec[4*i + 3] = obsVec[4*(i-1) + 3];
    }
    //std::cout << "leftProjection  = " << leftProjection << " rightProjection = " << rightProjection << std::endl;

    // Load the projected pixels into the output obseration vector
    obsVec[4*i + 0] = leftProjection [0]; // x
    obsVec[4*i + 1] = leftProjection [1]; // y
    obsVec[4*i + 2] = rightProjection[0];
    obsVec[4*i + 3] = rightProjection[1];
    
  } // End loop through points
  
  return obsVec;
}

}; // End class LrocPairModel
//===================================================================================================
//===================================================================================================

/// Functor to evaluate the residuals for left point observations (no camera parameters used)
struct LeftCostFunctor
{
private: // Variables
  LrocPairModel *_baseModel; // Make sure this does not go out of scope!
  double _observedRow;
  double _observedCol;
  
public:  // Functions
  
  /// Constructor
  LeftCostFunctor(LrocPairModel *baseModel, const double observedRow, const double observedCol)
  {
    _baseModel   = baseModel;
    _observedRow = observedRow;
    _observedCol = observedCol;
  }
  
  /// Wrapper for left observation function
  bool operator()(const double* const point, double* residuals) const 
  {
    double observations[2];
    if (!_baseModel->getLeftObservation(point, observations))
      return false;
    residuals[0] = observations[0] - _observedRow;
    residuals[1] = observations[1] - _observedCol;
    return true;
  }
}; // end struct LeftCostFunctor


/// Functor to evaluate the residuals for right point observations (includes camera parameters)
struct RightCostFunctor
{
private: // Variables
  LrocPairModel *_baseModel; // Make sure this does not go out of scope!
  double _observedRow;
  double _observedCol;
  
public:  // Functions
  
  /// Constructor
  RightCostFunctor(LrocPairModel *baseModel, const double observedRow, const double observedCol)
  {
    _baseModel   = baseModel;
    _observedRow = observedRow;
    _observedCol = observedCol;
  }
  
  /// Wrapper for left observation function
  bool operator()(const double* const camera, const double* const point, double* residuals) const 
  {
    double observations[2];
    if (!_baseModel->getRightObservation(camera, point, observations))
      return false;
    residuals[0] = observations[0] - _observedRow;
    residuals[1] = observations[1] - _observedCol;
    return true;
  }
}; // end struct RightCostFunctor







//----------------------------------------------------------------------------------------------------------------
