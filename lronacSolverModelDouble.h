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
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Stereo/PreFilter.h>
#include <vw/Stereo/Correlate.h>
#include <vw/Stereo/StereoModel.h>
#include <vw/Math/LevenbergMarquardt.h>
#include <vw/Math/EulerAngles.h>

#include <asp/IsisIO/IsisCameraModel.h> //::point_to_pixel>
#include <asp/Core/IntegralAutoGainDetector.h>
//#include <asp/IsisIO/IsisInterfaceLineScan.h>

#include <stereo.h>

#include <IsisInterfaceLineScanRot.h>
#include <lronacSolverSupport.h>


#include "ceres/ceres.h"
#include "glog/logging.h"

#ifndef __LRONACPIPELINE_SOLVERMODELDOUBLE_H__
#define __LRONACPIPELINE_SOLVERMODELDOUBLE_H__


// Points are passed in four sets of pairs.  Could try more overlap in the future.
struct PointObsList
{
  std::vector<vw::Vector2> leftObsList; // These need to be the same size
  std::vector<vw::Vector2> rightObsList;
  
  size_t size() const {return leftObsList.size();} ///< Return the number of points
};

/// Class for solving for the rotation between two LRONAC cameras
class LrocPairModel : public vw::math::LeastSquaresModelBase<LrocPairModel>
{
public:  // Definitions ------------------------------------------------------------------------------

/// * Defines a result_type that is the type returned by
///   evaluating the functor.  Typically Vector<float> or
///   Vector<double>
typedef vw::Vector<double> result_type;

/// * Defines a domain_type that is the type of the search
///   space.  Often a Vector<foo>, but can reflect other
///   topologies if needed.
typedef vw::Vector<double> domain_type;


/// * Defines a jacobian_type corresponding to the space of
///   jacobian matrices.  Typically Matrix<foo>.
typedef vw::Matrix<double> jacobian_type;

private: // Variables ---------------------------------------------------------------------------

  // Camera models
  IsisInterfaceLineScanRot _leftCameraModel;
  IsisInterfaceLineScanRot _rightCameraModel;

  IsisInterfaceLineScanRot _leftStereoCameraModel;
  IsisInterfaceLineScanRot _rightStereoCameraModel;

  mutable AdjustedCameraModelRot _leftStereoCameraRotatedModel;
  mutable AdjustedCameraModelRot _rightStereoCameraRotatedModel;

  // Observation records
  const PointObsList *_leftRight;   // Main pair
  const PointObsList *_leftSRightS; // Stereo pair
  const PointObsList *_leftLeftS;   // LE   to LE_S
  const PointObsList *_rightRightS; // RE   to RE_S
  const PointObsList *_leftRightS;  // LE   to RE_S
  const PointObsList *_leftSRight;  // LE_S to RE

  
public: // Functions  -----------------------------------------------------------------------------------

/// Constructor performs initialization
LrocPairModel(const std::string &leftCubePath,        const std::string &rightCubePath, 
              const std::string &leftStereoCubePath,  const std::string &rightStereoCubePath)
  : _leftCameraModel      (leftCubePath),       _rightCameraModel      (rightCubePath), 
    _leftStereoCameraModel(leftStereoCubePath), _rightStereoCameraModel(rightStereoCubePath),
    _leftStereoCameraRotatedModel (boost::shared_ptr<IsisInterfaceLineScanRot>(&_leftStereoCameraModel,
                                                                               boost::serialization::null_deleter()) ),
    _rightStereoCameraRotatedModel(boost::shared_ptr<IsisInterfaceLineScanRot>(&_rightStereoCameraModel,
                                                                               boost::serialization::null_deleter()) )
{
  // Both camera models are loaded from file on initialization
  printf("Done constructing LROC model\n");
}

size_t getNumPoints() const
{
  return _leftRight->size() + _leftSRightS->size() +
         _leftLeftS->size() + _rightRightS->size() +
         _leftRightS->size() + _leftSRight->size();
}

///// Test functions to get a single pixel vector
//Vector3 getLeftVector (const Vector2& pixel) {return vw::math::normalize(_leftCameraModel.pixel_to_vector (pixel));}
/*  Vector3 getRightVector(const Vector2& pixel) 
{
  if (_solveWorldFrame)
    return vw::math::normalize(_rightCameraRotatedModel.pixel_to_vector(pixel));
  else // Camera frame
    return vw::math::normalize(_rightCameraModel.pixel_to_vector(pixel));
}
*/
/*
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
*/



//left-right   (0-2 only)
//leftS-rightS (9-11 only)
//left-leftS   (3-8 only)
//right-rightS (ALL!)

/// Helper function used by getInitialStateEstimate
/// - Computes the best estimate of a point location given two observations
bool computePointLocation(const vw::camera::CameraModel &cam1,
                          const vw::camera::CameraModel &cam2,
                          const vw::Vector2 pixel1, const vw::Vector2 pixel2,
                          const bool useStereo,
                          vw::Vector3 &pointLocation)
{
  const double MOON_RADIUS_M    = 1737400.0;
  vw::Vector3 leftCamCenter  = cam1.camera_center(pixel1);
  vw::Vector3 rightCamCenter = cam2.camera_center(pixel2);
  
  //TODO: Need to incorporate local rotations?
  vw::Vector3 leftVec  = vw::math::normalize(cam1.pixel_to_vector(pixel1));
  vw::Vector3 rightVec = vw::math::normalize(cam2.pixel_to_vector(pixel2));

  vw::Vector3 v12 = cross_prod(leftVec, rightVec);
  vw::Vector3 v1  = cross_prod(v12,     leftVec);
  vw::Vector3 v2  = cross_prod(v12,     rightVec);

  vw::Vector3 pointGuess1, pointGuess2;
  
  if (useStereo)
  {
    pointGuess1 = leftCamCenter  + dot_prod(v2, rightCamCenter-leftCamCenter )/dot_prod(v2, leftVec )*leftVec;
    pointGuess2 = rightCamCenter + dot_prod(v1, leftCamCenter -rightCamCenter)/dot_prod(v1, rightVec)*rightVec;
  
    //if ((i % 100000) == 0) // Display progress
    //  printf("%d\n", i);
  }
  else // Use spheroid ground intersections instead
  {
    // Things don't work if this is done with an initial state
    if (!raySphereIntersect(leftCamCenter,  leftVec,  MOON_RADIUS_M, pointGuess1))
      printf("Failed to intersect moon with left!\n");
    if (!raySphereIntersect(rightCamCenter, rightVec, MOON_RADIUS_M, pointGuess2))
      printf("Failed to intersect moon with right!\n");
  }

  // Starting 3D location is the midpoint of where the two cameras place it
  pointLocation = 0.5 * (pointGuess1 + pointGuess2);

  // Code to assist in debugging
  vw::Vector3 errorVec = pointGuess1 - pointGuess2;
  double triangulationError = vw::math::norm_2(errorVec);
  //std::cout << "Pixel 1 = " << pixel1 << ", Pixel2 = " << pixel2 << ", triangulationError = " << triangulationError << std::endl;

  // Sanity check
  const double intersectionRadius = vw::math::norm_2(pointLocation);
  const double camRadius          = vw::math::norm_2(rightCamCenter);
  //std::cout << "intersection radius = " << intersectionRadius << std::endl;
  
  if (camRadius < intersectionRadius)
  {
    printf("Warning: Point reverse intersection!  Using previous point location as an estimate.\n");
    printf("%lf < %lf\n", camRadius, intersectionRadius);
    pointLocation = leftCamCenter + (leftVec*37000); // Extend 37km from left camera
  } 

  return true;
}

/// Given the pixel pair observations, compute the initial state estimate.
/// - This also loads the observation vectors and returns a packed version of them.
bool getInitialStateEstimate(const PointObsList &leftRight,   // Main pair
                             const PointObsList &leftSRightS, // Stereo pair
                             const PointObsList &leftLeftS,   // LE to LE
                             const PointObsList &rightRightS, // RE to RE
                             const PointObsList &leftRightS,  // LE to RE_S
                             const PointObsList &leftSRight,  // LE_S to RE
                             vw::Vector<double> &stateEstimate, //Vector<double> &packedObsVector,
                             const std::vector<double> &inputState = std::vector<double>())
{ 
  const size_t PARAMS_PER_POINT = 3;
  const double MOON_RADIUS_M    = 1737400.0;
  
  // Record observation vectors to class variables
  _leftRight   = &leftRight;   // Main   pair
  _leftSRightS = &leftSRightS; // Stereo pair
  _leftLeftS   = &leftLeftS;   // LE   to LE
  _rightRightS = &rightRightS; // RE   to RE
  _leftRightS  = &leftRightS;  // LE   to RE_S
  _leftSRight  = &leftSRight;  // LE_S to RE

  /* Input parameter set:
   0  rotationOffsetX (local rotations applied to RE camera)
   1  rotationOffsetY
   2  rotationOffsetZ
   3  rotationOffsetSX (planet-centered rotations applied to both stereo cameras)
   4  rotationOffsetSY
   5  rotationOffsetSZ
   6  positionOffsetX  (Position offset applied to the stereo pair)
   7  positionOffsetY
   8  positionOffsetZ
   9  rotationOffsetSRX (additional local rotations applied to stereo RE camera)
   10 rotationOffsetSRY
   11 rotationOffsetSRZ
      x                 [In planet coordinates, repeated for every correspondence point]
      y
      z
  */
  
  // Determine the number of state elements
  const size_t numMainPairs       = leftRight.size();
  const size_t numStereoPairs     = leftSRightS.size();
  const size_t numLeftPairs       = leftLeftS.size();
  const size_t numRightPairs      = rightRightS.size();
  const size_t numLeftCrossPairs  = leftRightS.size();
  const size_t numRightCrossPairs = leftSRight.size();
  const size_t numPoints          = this->getNumPoints();
  const size_t numCamParams       = 12;
  printf("numMainPairs       = %d\n", numMainPairs);
  printf("numStereoPairs     = %d\n", numStereoPairs);
  printf("numLeftPairs       = %d\n", numLeftPairs);
  printf("numRightPairs      = %d\n", numRightPairs);
  printf("numLeftCrossPairs  = %d\n", numLeftCrossPairs);
  printf("numRightCrossPairs = %d\n", numRightCrossPairs);
  printf("numPoints          = %d\n", numPoints);
  
  const size_t numStateElements = numPoints*PARAMS_PER_POINT + numCamParams;
  stateEstimate.set_size(numStateElements);
  //packedObsVector.set_size(numPoints*4);

  // Camera parameters start at zero
  for (size_t i=0; i<numCamParams; ++i)
    stateEstimate[i] = 0;
  
  // If an initial state was provided, use those values instead of zeroes.
  if (inputState.size() > 0)
  {
    printf("Calculating initial points using input state...\n");
    for (size_t i=0; i<inputState.size(); ++i)
      stateEstimate[i] = inputState[i];
  }

  // Set up initial state vectors
  vw::Vector3 globalRotVec(stateEstimate[3], stateEstimate[4], stateEstimate[5]);
  vw::Vector3 globalPosVec(stateEstimate[6], stateEstimate[7], stateEstimate[8]);
  
  // Apply global transformation to the stereo pair
  _leftStereoCameraRotatedModel.set_axis_angle_rotation (globalRotVec);
  _rightStereoCameraRotatedModel.set_axis_angle_rotation(globalRotVec);
  _leftStereoCameraRotatedModel.set_translation (globalPosVec);
  _rightStereoCameraRotatedModel.set_translation(globalPosVec);

  //DEBUG
  // Set up georeference class with default moon datum
  vw::cartography::Datum datum("D_MOON");

  vw::cartography::GeoReference leftGeoRef(datum);
  vw::cartography::GeoReference rightGeoRef(datum);

  printf("Setting up state estimate\n");

  const bool useStereo = (inputState.size() > 0);
  vw::Vector3 pointLoc;
  
  // LEFT-RIGHT comparison
  size_t pointOffset = numCamParams; // Where to start writing the points in the state estimate
  for (size_t i=0; i<numMainPairs; ++i)
  {
    vw::Vector2 leftPixel  = leftRight.leftObsList[i];
    vw::Vector2 rightPixel = leftRight.rightObsList[i];
    computePointLocation(_leftCameraModel, _rightCameraModel, leftPixel, rightPixel, useStereo, pointLoc);
    
    // Record the x/y/z value for this point
    for (size_t p=0; p<PARAMS_PER_POINT; ++p)
      stateEstimate[pointOffset + i*PARAMS_PER_POINT + p] = pointLoc[p];
  } // End loop through points
  pointOffset += numMainPairs*PARAMS_PER_POINT;
  
  // LEFT STEREO-RIGHT STEREO comparison
  for (size_t i=0; i<numStereoPairs; ++i)
  {
    vw::Vector2 leftPixel  = leftSRightS.leftObsList[i];
    vw::Vector2 rightPixel = leftSRightS.rightObsList[i];
    computePointLocation(_leftCameraModel, _rightCameraModel, leftPixel, rightPixel, useStereo, pointLoc);
    
    // Record the x/y/z value for this point
    for (size_t p=0; p<PARAMS_PER_POINT; ++p)
      stateEstimate[pointOffset + i*PARAMS_PER_POINT + p] = pointLoc[p];
  } // End loop through points
  pointOffset += numStereoPairs*PARAMS_PER_POINT;
  
  // LEFT-LEFT STEREO comparison
  for (size_t i=0; i<numLeftPairs; ++i)
  {
    vw::Vector2 leftPixel  = leftLeftS.leftObsList[i];
    vw::Vector2 rightPixel = leftLeftS.rightObsList[i];
    computePointLocation(_leftCameraModel, _leftStereoCameraRotatedModel, leftPixel, rightPixel, useStereo, pointLoc);
    
    // Record the x/y/z value for this point
    for (size_t p=0; p<PARAMS_PER_POINT; ++p)
      stateEstimate[pointOffset + i*PARAMS_PER_POINT + p] = pointLoc[p];
  } // End loop through points
  pointOffset += numLeftPairs*PARAMS_PER_POINT;
  
  // RIGHT-RIGHT STEREO comparison
  for (size_t i=0; i<numRightPairs; ++i)
  {
    vw::Vector2 leftPixel  = rightRightS.leftObsList[i];
    vw::Vector2 rightPixel = rightRightS.rightObsList[i];
    computePointLocation(_rightCameraModel, _rightStereoCameraRotatedModel, leftPixel, rightPixel, useStereo, pointLoc);
    
    // Record the x/y/z value for this point
    for (size_t p=0; p<PARAMS_PER_POINT; ++p)
      stateEstimate[pointOffset + i*PARAMS_PER_POINT + p] = pointLoc[p];
  } // End loop through points
  pointOffset += numRightPairs*PARAMS_PER_POINT;
  
  // Now handle the two optional point sections

  if (numLeftCrossPairs > 0)
  {
    // LEFT-RIGHTS comparison
    for (size_t i=0; i<numLeftCrossPairs; ++i)
    {
      vw::Vector2 leftPixel  = leftRightS.leftObsList[i];
      vw::Vector2 rightPixel = leftRightS.rightObsList[i];
      computePointLocation(_leftCameraModel, _rightStereoCameraRotatedModel, leftPixel, rightPixel, useStereo, pointLoc);

      // Record the x/y/z value for this point
      for (size_t p=0; p<PARAMS_PER_POINT; ++p)
        stateEstimate[pointOffset + i*PARAMS_PER_POINT + p] = pointLoc[p];
    } // End loop through points
    pointOffset += numLeftCrossPairs*PARAMS_PER_POINT;
  }

  if (numRightCrossPairs > 0)
  {
    // LEFTS-RIGHT comparison
    for (size_t i=0; i<numRightCrossPairs; ++i)
    {
      vw::Vector2 leftPixel  = leftSRight.leftObsList[i];
      vw::Vector2 rightPixel = leftSRight.rightObsList[i];
      computePointLocation(_leftStereoCameraRotatedModel, _rightCameraModel, leftPixel, rightPixel, useStereo, pointLoc);

      // Record the x/y/z value for this point
      for (size_t p=0; p<PARAMS_PER_POINT; ++p)
        stateEstimate[pointOffset + i*PARAMS_PER_POINT + p] = pointLoc[p];
    } // End loop through points
    pointOffset += numRightCrossPairs*PARAMS_PER_POINT;
  }


  // At this point the state estimate is fully populated
  
  
  return true;
} // end getInitialStateEstimate()


/// Generates a vector containing the stereo intersection distance for each pixel pair
/// - This gives a good indication of how accurate the current correction is
std::vector<double> computeError(domain_type const& x)
{    
  // Determine the number of state elements
  const size_t numPoints = this->getNumPoints();
  std::vector<double> errorVector(numPoints);

  // Determine the number of state elements
  const size_t numMainPairs       = _leftRight->size();
  const size_t numStereoPairs     = _leftSRightS->size();
  const size_t numLeftPairs       = _leftLeftS->size();
  const size_t numRightPairs      = _rightRightS->size();
  const size_t numLeftCrossPairs  = _leftRightS->size();
  const size_t numRightCrossPairs = _leftSRight->size();
  const size_t numCamParams       = 12;

  // Set up initial state vectors
  vw::Vector3 localRotVec      (x[0], x[1], x[2]);
  vw::Vector3 globalRotVec     (x[3], x[4], x[5]);
  vw::Vector3 globalPosVec     (x[6], x[7], x[8]);
  vw::Vector3 stereoLocalRotVec(x[9], x[10], x[11]);
  vw::Vector3 nullVec(0,0,0); // Just to use our custom rotation function

  // Apply global transformation to the stereo pair
  _leftStereoCameraRotatedModel.set_axis_angle_rotation(globalRotVec);
  _rightStereoCameraRotatedModel.set_axis_angle_rotation(globalRotVec);
  _leftStereoCameraRotatedModel.set_translation(globalPosVec);
  _rightStereoCameraRotatedModel.set_translation(globalPosVec);

  const int PARAMS_PER_POINT = 3;

  // LEFT-RIGHT comparison
  size_t pointOffset = numCamParams; // Where to start writing the points in the state estimate
  size_t errorOffset = 0;
  for (size_t i=0; i<numMainPairs; ++i)
  {
    // Create a point object
    vw::Vector3 thisPoint(x[pointOffset + i*PARAMS_PER_POINT + 0],
                          x[pointOffset + i*PARAMS_PER_POINT + 1],
                          x[pointOffset + i*PARAMS_PER_POINT + 2]);

    // Project to pixel locations
    vw::Vector2 leftProjection  = _leftCameraModel.point_to_pixel_rotated (thisPoint, nullVec,     _leftRight->leftObsList[i][1]);
    vw::Vector2 rightProjection = _rightCameraModel.point_to_pixel_rotated(thisPoint, localRotVec, _leftRight->rightObsList[i][1]); // This treats the angles as euler angles!

    // Compare to observed pixel locations
    vw::Vector2 leftDiff      = leftProjection  - _leftRight->leftObsList [i];
    vw::Vector2 rightDiff     = rightProjection - _leftRight->rightObsList[i];
    double  leftError     = vw::math::norm_2(leftDiff );
    double  rightError    = vw::math::norm_2(rightDiff);
    errorVector[errorOffset + i] = (leftError + rightError) / 2.0; // Returned value is average of left and right position error
  } // End loop through points
  pointOffset += numMainPairs*PARAMS_PER_POINT;
  errorOffset += numMainPairs;

  // LEFT STEREO-RIGHT STEREO comparison
  for (size_t i=0; i<numStereoPairs; ++i)
  {
    // Create a point object
    vw::Vector3 thisPoint(x[pointOffset + i*PARAMS_PER_POINT + 0],
                      x[pointOffset + i*PARAMS_PER_POINT + 1],
                      x[pointOffset + i*PARAMS_PER_POINT + 2]);

    // Project to pixel locations
    vw::Vector2 leftProjection  = _leftStereoCameraRotatedModel.point_to_pixel_rotated (thisPoint, nullVec,           _leftSRightS->leftObsList[i][1]);
    vw::Vector2 rightProjection = _rightStereoCameraRotatedModel.point_to_pixel_rotated(thisPoint, stereoLocalRotVec, _leftSRightS->rightObsList[i][1]); // This treats the angles as euler angles!

    // Compare to observed pixel locations
    vw::Vector2 leftDiff      = leftProjection  - _leftSRightS->leftObsList [i];
    vw::Vector2 rightDiff     = rightProjection - _leftSRightS->rightObsList[i];
    double  leftError     = vw::math::norm_2(leftDiff );
    double  rightError    = vw::math::norm_2(rightDiff);
    errorVector[errorOffset + i] = (leftError + rightError) / 2.0; // Returned value is average of left and right position error
  } // End loop through points
  pointOffset += numStereoPairs*PARAMS_PER_POINT;
  errorOffset += numStereoPairs;

  // LEFT-LEFT STEREO comparison
  for (size_t i=0; i<numLeftPairs; ++i)
  {
    // Create a point object
    vw::Vector3 thisPoint(x[pointOffset + i*PARAMS_PER_POINT + 0],
                      x[pointOffset + i*PARAMS_PER_POINT + 1],
                      x[pointOffset + i*PARAMS_PER_POINT + 2]);

    // Project to pixel locations
    vw::Vector2 leftProjection  = _leftCameraModel.point_to_pixel_rotated             (thisPoint, nullVec, _leftLeftS->leftObsList[i][1]);
    vw::Vector2 rightProjection = _leftStereoCameraRotatedModel.point_to_pixel_rotated(thisPoint, nullVec, _leftLeftS->rightObsList[i][1]); // This treats the angles as euler angles!

    //std::cout << "Point: " << thisPoint << " projected to " << leftProjection << ", " << rightProjection << std::endl;

    // Compare to observed pixel locations
    vw::Vector2 leftDiff      = leftProjection  - _leftLeftS->leftObsList [i];
    vw::Vector2 rightDiff     = rightProjection - _leftLeftS->rightObsList[i];
    double  leftError     = vw::math::norm_2(leftDiff );
    double  rightError    = vw::math::norm_2(rightDiff);
    errorVector[errorOffset + i] = (leftError + rightError) / 2.0; // Returned value is average of left and right position error
  } // End loop through points
  pointOffset += numLeftPairs*PARAMS_PER_POINT;
  errorOffset += numLeftPairs;

  //std::cout << "rotVec = " << globalRotVec << ", offsetVec = " << globalPosVec << std::endl;

  // RIGHT-RIGHT STEREO comparison
  for (size_t i=0; i<numRightPairs; ++i)
  {
    // Create a point object
    vw::Vector3 thisPoint(x[pointOffset + i*PARAMS_PER_POINT + 0],
                      x[pointOffset + i*PARAMS_PER_POINT + 1],
                      x[pointOffset + i*PARAMS_PER_POINT + 2]);

    // Project to pixel locations
    vw::Vector2 leftProjection  = _rightCameraModel.point_to_pixel_rotated             (thisPoint, localRotVec,       _rightRightS->leftObsList[i][1]);
    vw::Vector2 rightProjection = _rightStereoCameraRotatedModel.point_to_pixel_rotated(thisPoint, stereoLocalRotVec, _rightRightS->rightObsList[i][1]); // This treats the angles as euler angles!

    // Compare to observed pixel locations
    vw::Vector2 leftDiff      = leftProjection  - _rightRightS->leftObsList [i];
    vw::Vector2 rightDiff     = rightProjection - _rightRightS->rightObsList[i];
    double  leftError     = vw::math::norm_2(leftDiff );
    double  rightError    = vw::math::norm_2(rightDiff);
    errorVector[errorOffset + i] = (leftError + rightError) / 2.0; // Returned value is average of left and right position error
  } // End loop through points
  pointOffset += numRightPairs*PARAMS_PER_POINT;
  errorOffset += numRightPairs;

  // Now handle the optional point pairs

  if (numLeftCrossPairs > 0)
  {
    // LEFT-RIGHT STEREO
    for (size_t i=0; i<numLeftCrossPairs; ++i)
    {
      // Create a point object
      vw::Vector3 thisPoint(x[pointOffset + i*PARAMS_PER_POINT + 0],
                        x[pointOffset + i*PARAMS_PER_POINT + 1],
                        x[pointOffset + i*PARAMS_PER_POINT + 2]);

      // Project to pixel locations
      vw::Vector2 leftProjection  = _leftCameraModel.point_to_pixel_rotated              (thisPoint, nullVec,           _leftRightS->leftObsList[i][1]);
      vw::Vector2 rightProjection = _rightStereoCameraRotatedModel.point_to_pixel_rotated(thisPoint, stereoLocalRotVec, _leftRightS->rightObsList[i][1]); // This treats the angles as euler angles!

      // Compare to observed pixel locations
      vw::Vector2 leftDiff      = leftProjection  - _leftRightS->leftObsList [i];
      vw::Vector2 rightDiff     = rightProjection - _leftRightS->rightObsList[i];
      double  leftError     = vw::math::norm_2(leftDiff );
      double  rightError    = vw::math::norm_2(rightDiff);
      errorVector[errorOffset + i] = (leftError + rightError) / 2.0; // Returned value is average of left and right position error
    } // End loop through points
    pointOffset += numLeftCrossPairs*PARAMS_PER_POINT;
    errorOffset += numLeftCrossPairs;
  }

  if (numRightCrossPairs > 0)
  {
    // LEFT STEREO-RIGHT
    for (size_t i=0; i<numRightCrossPairs; ++i)
    {
      // Create a point object
      vw::Vector3 thisPoint(x[pointOffset + i*PARAMS_PER_POINT + 0],
                        x[pointOffset + i*PARAMS_PER_POINT + 1],
                        x[pointOffset + i*PARAMS_PER_POINT + 2]);

      // Project to pixel locations
      vw::Vector2 leftProjection  = _leftStereoCameraRotatedModel.point_to_pixel_rotated(thisPoint, nullVec,     _leftSRight->leftObsList[i][1]);
      vw::Vector2 rightProjection = _rightCameraModel.point_to_pixel_rotated            (thisPoint, localRotVec, _leftSRight->rightObsList[i][1]); // This treats the angles as euler angles!

      // Compare to observed pixel locations
      vw::Vector2 leftDiff      = leftProjection  - _leftSRight->leftObsList [i];
      vw::Vector2 rightDiff     = rightProjection - _leftSRight->rightObsList[i];
      double  leftError     = vw::math::norm_2(leftDiff );
      double  rightError    = vw::math::norm_2(rightDiff);
      errorVector[errorOffset + i] = (leftError + rightError) / 2.0; // Returned value is average of left and right position error
    } // End loop through points
    pointOffset += numRightCrossPairs*PARAMS_PER_POINT;
    errorOffset += numRightCrossPairs;
  }

  return errorVector;
    
} // end computeError()

/// For a given point, compute the left camera observation.
/// - Returns false if the point is not visible.
bool getLeftObservation(const double* const pointParams, double *observation, int guessRow=-1) const
{
  // Create a point object
  vw::Vector3 thisPoint(pointParams[0], pointParams[1], pointParams[2]);
  vw::Vector2 projection;
  try // Project the point into the camera
  {
    vw::Vector3 nullVec(0,0,0); // Left camera not currently rotated
    projection  = _leftCameraModel.point_to_pixel_rotated(thisPoint, nullVec, guessRow);
  }
  catch(std::exception& e) // Handle errors
  {
    std::cout << "Warning: Failed to project location: " << thisPoint << ", what= ." << e.what() << std::endl;
    return false;
  }
  
  // Load the projected pixels into the output observation vector
  observation[0] = projection[0]; // x
  observation[1] = projection[1]; // y

  return true;
}

/// For a given point, compute the right camera observation.
/// - Returns false if the point is not visible.
/// - rotAngleParams points to parameters 0-2
bool getRightObservation(const double* const rotAngleParams, const double* const pointParams, double *observation, int guessRow=-1)
{
  // This function returns an error vector for a given set of parameters

  // Apply the rotations from the state vector to the right LROC camera model
  vw::Vector3 localRotVec(rotAngleParams[0], rotAngleParams[1], rotAngleParams[2]);

  // Create a point object
  vw::Vector3 thisPoint(pointParams[0], pointParams[1], pointParams[2]);
  vw::Vector2 projection;
  try // Project the point into the camera
  {
    projection = _rightCameraModel.point_to_pixel_rotated(thisPoint, localRotVec, guessRow);
  }
  catch(std::exception& e)
  {
    std::cout << "Warning: Failed to project location: " << thisPoint << ", what= ." << e.what() << std::endl;
    return false;
  }

  // Load the projected pixels into the output obserVation vector
  observation[0] = projection[0];
  observation[1] = projection[1];

  return true;
}

/// For a given point, compute the left Stereo camera observation.
/// - Returns false if the point is not visible.
/// - rotParams points to parameters 3-5
/// - posParams points to parameters 6-8
bool getLeftStereoObservation(const double* const rotParams, const double* const posParams, const double* const pointParams, double *observation, int guessRow=-1)
{
  // This function returns an error vector for a given set of parameters

  // Apply the rotations from the state vector to the right LROC camera model
  vw::Vector3 rotVec   (rotParams[0], rotParams[1], rotParams[2]);
  vw::Vector3 offsetVec(posParams[0], posParams[1], posParams[2]);
  _leftStereoCameraRotatedModel.set_axis_angle_rotation(rotVec);
  _leftStereoCameraRotatedModel.set_translation(offsetVec);

  // Create a point object
  vw::Vector3 thisPoint(pointParams[0], pointParams[1], pointParams[2]);
  vw::Vector2 projection;
  try // Project the point into the camera
  {
    vw::Vector3 nullVec(0,0,0); // Left camera not currently rotated
    projection = _leftStereoCameraRotatedModel.point_to_pixel_rotated(thisPoint, nullVec, guessRow);
  }
  catch(std::exception& e)
  {
    std::cout << "Warning: Failed to project location: " << thisPoint << ", what= ." << e.what() << std::endl;
    return false;
  }

  // Load the projected pixels into the output obserVation vector
  observation[0] = projection[0];
  observation[1] = projection[1];

  return true;
}

/// For a given point, compute the right Stereo camera observation.
/// - Returns false if the point is not visible.
/// - rotParams points to parameters 3-5
/// - posParams points to parameters 6-8
/// - localRotParams points to parameters 9-11
bool getRightStereoObservation(const double* const rotParams, 
                               const double* const posParams,
                               const double* const localRotParams,
                               const double* const pointParams, double *observation, int guessRow=-1)
{
  // This function returns an error vector for a given set of parameters

  // Apply the rotations from the state vector to the right LROC camera model
  vw::Vector3 rotVec     (rotParams[0],      rotParams[1],      rotParams[2]);
  vw::Vector3 offsetVec  (posParams[0],      posParams[1],      posParams[2]);
  vw::Vector3 localRotVec(localRotParams[0], localRotParams[1], localRotParams[2]);
  _rightStereoCameraRotatedModel.set_axis_angle_rotation(rotVec);
  _rightStereoCameraRotatedModel.set_translation(offsetVec);

  // Create a point object
  vw::Vector3 thisPoint(pointParams[0], pointParams[1], pointParams[2]);
  vw::Vector2 projection;
  try // Project the point into the camera
  {
    projection = _rightStereoCameraRotatedModel.point_to_pixel_rotated(thisPoint, localRotVec, guessRow);
  }
  catch(std::exception& e)
  {
    std::cout << "Warning: Failed to project location: " << thisPoint << ", what= ." << e.what() << std::endl;
    return false;
  }

  // Load the projected pixels into the output obserVation vector
  observation[0] = projection[0];
  observation[1] = projection[1];

  return true;
}



/// * Defines a method: result_type operator()( domain_type const& x ) const;
///   that evaluates the model function at the given point.
result_type operator()( domain_type const& x ) const
{
  // This function returns an error vector for a given set of parameters

  // Apply the rotations from the state vector to the right LROC camera model
  vw::Vector3 rotVec   (x[0], x[1], x[2]);
  vw::Vector3 offsetVec(x[3], x[4], x[5]); // Only used if _includePosition = true

  //TODO: This stuff is out of date
  _rightStereoCameraRotatedModel.set_translation(offsetVec);
  _rightStereoCameraRotatedModel.set_axis_angle_rotation(rotVec);
  
  const double rad2deg = 180.0 / M_PI;
  //printf("Trying rotation %lf, %lf, %lf\n", x[0]*rad2deg, x[1]*rad2deg, x[2]*rad2deg);
  
  //std::cout << "Left  camera center = " << _leftCameraModel.camera_center()  << std::endl;
  //std::cout << "Right camera center = " << _rightCameraModel.camera_center() << std::endl;
  
  // Set up output vector
  const size_t numPoints = this->getNumPoints();
  vw::Vector<double> obsVec(numPoints*4);
  
  //TODO: Probably don't need to update this function since we are using CERES instead
  
  /*
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
      Vector3 nullVec(0,0,0);
      leftProjection  = _leftCameraModel.point_to_pixel_rotated(thisPoint, nullVec, _leftRows[i]);
      if (_solveWorldFrame)
        rightProjection = _rightCameraRotatedModel.point_to_pixel(thisPoint);
      else // Camera frame
        rightProjection = _rightCameraModel.point_to_pixel_rotated(thisPoint, rotVec, _rightRows[i]);
    }
    catch(...)
    {
      std::cout << "Warning: Failed to project location " << i << ": " << thisPoint << ", using previous intersection location X_X." << std::endl;
      
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
  */
  
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
  vw::Vector2 _observation;
  
public:  // Functions
  
  /// Constructor
  LeftCostFunctor(LrocPairModel *baseModel, const vw::Vector2 observation)
  {
    _baseModel   = baseModel;
    _observation = observation;
  }
  
  /// Wrapper for left observation function
  bool operator()(const double* const point, double* residuals) const 
  {
    double observations[2];
    if (!_baseModel->getLeftObservation(point, observations, _observation[1]))
      return false;
    //std::ofstream file("/home/smcmich1/logDoubleL.csv", std::ofstream::app);
    ////printf("observations = %lf, %lf\n", observations[0], observations[1]);
    //file << observations[0] << ", " << observations[1] << std::endl;
    //file.close();
    residuals[0] = observations[0] - _observation[0];
    residuals[1] = observations[1] - _observation[1];
    return true;
  }
}; // end struct LeftCostFunctor


/// Functor to evaluate the residuals for right point observations (includes camera parameters)
struct RightCostFunctor
{
private: // Variables
  LrocPairModel *_baseModel; // Make sure this does not go out of scope!
  vw::Vector2 _observation;
  
public:  // Functions
  
  /// Constructor
  RightCostFunctor(LrocPairModel *baseModel, const vw::Vector2 observation)
  {
    _baseModel   = baseModel;
    _observation = observation;
  }
  
  /// Wrapper for right observation function
  bool operator()(const double* const camera, const double* const point, double* residuals) const 
  {
    double observations[2];
    if (!_baseModel->getRightObservation(camera, point, observations, _observation[1]))
      return false;
    residuals[0] = observations[0] - _observation[0];
    residuals[1] = observations[1] - _observation[1];
    return true;
  }
}; // end struct RightCostFunctor

/// Functor to evaluate the residuals for left stereo point observations (includes camera parameters)
struct LeftStereoCostFunctor
{
private: // Variables
  LrocPairModel *_baseModel; // Make sure this does not go out of scope!
  vw::Vector2 _observation;
  
public:  // Functions
  
  /// Constructor
  LeftStereoCostFunctor(LrocPairModel *baseModel, const vw::Vector2 observation)
  {
    _baseModel   = baseModel;
    _observation = observation;
  }
  
  /// Wrapper for left stereo observation function
  bool operator()(const double* const rotParams, const double* const posParams, const double* const point, double* residuals) const
  //bool operator()(const double* const rotParams, const double* const point, double* residuals) const
  {
    //const double* const posParams = &(rotParams[3]);
    double observations[2];
    if (!_baseModel->getLeftStereoObservation(rotParams, posParams, point, observations, _observation[1]))
      return false;
    //std::ofstream file("/home/smcmich1/logDoubleLS.csv", std::ofstream::app);
    ////printf("observations = %lf, %lf\n", observations[0], observations[1]);
    //file << observations[0] << ", " << observations[1] << std::endl;
    //file.close();
    residuals[0] = observations[0] - _observation[0];
    residuals[1] = observations[1] - _observation[1];
    return true;
  }
}; // end struct LeftStereoCostFunctor


/// Functor to evaluate the residuals for right stereo point observations (includes camera parameters)
struct RightStereoCostFunctor
{
private: // Variables
  LrocPairModel *_baseModel; // Make sure this does not go out of scope!
  vw::Vector2 _observation;
  
public:  // Functions
  
  /// Constructor
  RightStereoCostFunctor(LrocPairModel *baseModel, const vw::Vector2 observation)
  {
    _baseModel   = baseModel;
    _observation = observation;
  }
  
  /// Wrapper for right stereo observation function
  bool operator()(const double* const rotParams, 
                  const double* const posParams, 
                  const double* const localRotParams, 
                  const double* const point, double* residuals) const 
  {
    double observations[2];
    if (!_baseModel->getRightStereoObservation(rotParams, posParams, localRotParams, point, observations, _observation[1]))
      return false;
    residuals[0] = observations[0] - _observation[0];
    residuals[1] = observations[1] - _observation[1];
    return true;
  }
}; // end struct RightStereoCostFunctor


#endif


//----------------------------------------------------------------------------------------------------------------
