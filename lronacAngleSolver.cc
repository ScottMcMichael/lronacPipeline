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


/// \file lronacAngleSolver.cc
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

//#include <asp/Tools/stereo.h>
#include <stereo.h> // Copied to local directory


using namespace vw;
using namespace vw::stereo;
using namespace asp::isis;
using namespace asp;
using std::endl;
using std::setprecision;
using std::setw;




/*
//-----------------------------------------------------------------------------------------------------------
// Copy of class to allow easy insertion of debug statements
  template <class ImplT>
  typename ImplT::domain_type levenberg_marquardt_V( LeastSquaresModelBase<ImplT> const& least_squares_model,
                                                   typename ImplT::domain_type const& seed,
                                                   typename ImplT::result_type const& observation,
                                                   int &status,
                                                   double abs_tolerance = VW_MATH_LM_ABS_TOL,
                                                   double rel_tolerance = VW_MATH_LM_REL_TOL,
                                                   double max_iterations = 100) {

    status = optimization::eDidNotConverge;

    const ImplT& model = least_squares_model.impl();
    bool done = false;
    double Rinv = 10;
    double lambda = 0.1;

    typename ImplT::domain_type x_try, x = seed;
    typename ImplT::result_type h = model(x);
    typename ImplT::result_type error = model.difference(observation, h);
    double norm_start = norm_2(error);

    VW_OUT(DebugMessage, "math") << "LM: initial guess for the model is " << seed << std::endl;
    VW_OUT(VerboseDebugMessage, "math") << "LM: starting error " << error << std::endl;
    std::cout << "LM: starting error " << error << std::endl;
    VW_OUT(DebugMessage, "math") << "LM: starting norm is: " << norm_start << std::endl;

    // Solution may already be good enough
    if (norm_start < abs_tolerance) {
      status = optimization::eConvergedAbsTolerance;
      VW_OUT(DebugMessage, "math") << "CONVERGED TO ABSOLUTE TOLERANCE\n";
      done = true;
    }

    int outer_iter = 0;
    while (!done){

      bool shortCircuit = false;
      outer_iter++;
      VW_OUT(DebugMessage, "math") << "LM: outer iteration " << outer_iter << "   x = " << x << std::endl;
      std::cout << "LM: outer iteration " << outer_iter << "   x = " << x << std::endl;

      // Compute the value, derivative, and hessian of the cost function
      // at the current point.  These remain valid until the parameter
      // vector changes.

      // expected measurement with new x
      h = model(x);

      // Difference between observed and predicted and error (2-norm of difference)
      error = model.difference(observation, h);
      norm_start = norm_2(error);
      //VW_OUT(DebugMessage, "math") << "LM: outer iteration starting robust norm: " << norm_start << std::endl;

      // Measurement Jacobian
      typename ImplT::jacobian_type J = model.jacobian(x);

      Vector<double> del_J = -1.0 * Rinv * (transpose(J) * error);

      std::ofstream delJFile("/home/smcmich1/del_J.csv"); 
      for (size_t i=0; i<del_J.size(); ++i)
        delJFile << del_J[i] << std::endl;
      delJFile.close();
      
      // Hessian of cost function (using Gauss-Newton approximation)
      Matrix<double> hessian = Rinv * (transpose(J) * J);

      int iterations = 0;
      double norm_try = norm_start+1.0;
      while (norm_try > norm_start){

        // Increase diagonal elements to dynamically mix gradient
        // descent and Gauss-Newton.
        Matrix<double> hessian_lm = hessian;
        for ( unsigned i=0; i < hessian_lm.rows(); ++i ){
          hessian_lm(i,i) += hessian_lm(i,i)*lambda + lambda;
        }

        // Solve for update
        typename ImplT::domain_type delta_x;
        if (hessian_lm.rows() <= 2 && det(hessian_lm) > 0.0){
          // Direct method is more efficient for small matrices, also
          // here we avoid calling LAPACK which we've seen misbehave
          // in this situation in a multi-threaded environment.
          delta_x = inverse(hessian_lm)*del_J;
        }else{
          try{
            // By construction, hessian_lm is symmetric and
            // positive-definite.
            delta_x = solve_symmetric(hessian_lm, del_J);
          }catch ( const ArgumentErr& e ) {
            // If lambda is very small, the matrix becomes numerically
            // singular. In that case use the more general
            // least_squares solver.
            delta_x = least_squares(hessian_lm, del_J);
          }

        }

        // update parameter vector
        x_try = x - delta_x;

        typename ImplT::result_type h_try = model(x_try);

        typename ImplT::result_type error_try = model.difference(observation, h_try);
        norm_try = norm_2(error_try);

        //VW_OUT(VerboseDebugMessage, "math") << "LM: inner iteration " << iterations << " error is " << error_try << std::endl;
        //VW_OUT(DebugMessage, "math") << "\tLM: inner iteration " << iterations << " norm is " << norm_try << std::endl;
        std::cout << "\tLM: inner iteration " << iterations << " norm is " << norm_try << std::endl;

        if (norm_try > norm_start)
          // Increase lambda and try again
          lambda *= 10;

        ++iterations; // Sanity check on iterations in this loop
        if (iterations > 20) {
          //VW_OUT(DebugMessage, "math") << "\n****LM: too many inner iterations - short circuiting\n" << std::endl;
          std::cout << "\n****LM: too many inner iterations - short circuiting\n" << std::endl;
          shortCircuit = true;
          norm_try = norm_start;
        }
        //VW_OUT(DebugMessage, "math") << "\tlambda = " << lambda << std::endl;
        std::cout << "\tlambda = " << lambda << std::endl;
        
        std::ofstream deltaXFile("/home/smcmich1/deltaX.csv"); 
        for (size_t i=0; i<delta_x.size(); ++i)
          deltaXFile << delta_x[i] << std::endl;
        deltaXFile.close();

        std::ofstream finalHessianFile("/home/smcmich1/hessian.csv");
        finalHessianFile.precision(3);
        for (int r=0; r<hessian_lm.rows(); ++r)
        {
          for (int c=0; c<hessian_lm.cols(); ++c)
          {
            finalHessianFile << hessian_lm(r,c) << ", ";
          }
          finalHessianFile << std::endl;
        }
        finalHessianFile.close();
        
      } // End inner loop

      // Percentage change convergence criterion
      if (((norm_start-norm_try)/norm_start) < rel_tolerance) {
        status = optimization::eConvergedRelTolerance;
        VW_OUT(DebugMessage, "math") << "CONVERGED TO RELATIVE TOLERANCE\n";
        std::cout << "CONVERGED TO RELATIVE TOLERANCE\n";
        done = true;
      }

      // Absolute error convergence criterion
      if (norm_try < abs_tolerance) {
        status = optimization::eConvergedAbsTolerance;
        VW_OUT(DebugMessage, "math") << "CONVERGED TO ABSOLUTE TOLERANCE\n";
        std::cout << "CONVERGED TO ABSOLUTE TOLERANCE\n";
        done = true;
      }

      // Max iterations convergence criterion
      if (outer_iter >= max_iterations) {
        VW_OUT(DebugMessage, "math") << "REACHED MAX ITERATIONS!";
        std::cout << "REACHED MAX ITERATIONS!";
        done = true;
      }

      // Take trial parameters as new parameters
      // If we short-circuited the inner loop, then we didn't actually find a
      // better p, so don't update it.
      if (!shortCircuit)
        x = x_try;

      // Take trial error as new error
      norm_start = norm_try;

      // Decrease lambda
      lambda /= 10;
      //VW_OUT(DebugMessage, "math") << "lambda = " << lambda << std::endl;
      //VW_OUT(DebugMessage, "math") << "LM: end of outer iteration " << outer_iter << " with error " << norm_try << std::endl;
    
    } // End outer loop
    
    
    VW_OUT(DebugMessage, "math") << "LM: finished with: " << outer_iter << "\n";
    std::cout << "LM: finished with: " << outer_iter << "\n";
    return x;
  }

}} // end vw::math
*/

namespace asp {
namespace isis {

class IsisInterfaceLineScanRot : public asp::isis::IsisInterfaceLineScan 
{
public:
  IsisInterfaceLineScanRot( std::string const& file ) : asp::isis::IsisInterfaceLineScan(file) {}
  virtual ~IsisInterfaceLineScanRot() {}

  /// Additional function to apply an in-camera rotation during this operation
  vw::Vector2
    point_to_pixel_rotated( vw::Vector3 const& point, vw::Vector3 const& rotAngles) const;

};


/// Helper class for computing the projection line - With added rotation argument!
class EphemerisLMA_rot : public vw::math::LeastSquaresModelBase<EphemerisLMA_rot> {
  vw::Vector3 m_point;
  Isis::Camera* m_camera;
  Isis::CameraDistortionMap *m_distortmap;
  Isis::CameraFocalPlaneMap *m_focalmap;
  vw::Vector3 m_rotation;
public:
  typedef vw::Vector<double> result_type; // Back project result
  typedef vw::Vector<double> domain_type; // Ephemeris time
  typedef vw::Matrix<double> jacobian_type;

  inline EphemerisLMA_rot( vw::Vector3 const& point, vw::Vector3 const& rotation,
                       Isis::Camera* camera,
                       Isis::CameraDistortionMap* distortmap,
                       Isis::CameraFocalPlaneMap* focalmap ) 
                        : m_point(point), m_camera(camera), m_distortmap(distortmap), 
                          m_focalmap(focalmap), m_rotation(rotation) 
                       {}
  inline result_type operator()( domain_type const& x ) const;
};


// LMA for projecting point to linescan camera
EphemerisLMA_rot::result_type
EphemerisLMA_rot::operator()( EphemerisLMA_rot::domain_type const& x ) const 
{

  // Setting Ephemeris Time
  m_camera->setTime(Isis::iTime(x[0]));

  // Calculating the look direction in camera frame
  Vector3 instru;
  m_camera->instrumentPosition(&instru[0]);
  instru *= 1000;  // Spice gives in km, convert to meters
  
  Vector3 lookB = normalize( m_point - instru ); // Point to camera vector in body coordinate frame
  std::vector<double> lookB_copy(3);
  std::copy( lookB.begin(), lookB.end(), lookB_copy.begin() );
  
  std::vector<double> lookJ = m_camera->bodyRotation()->J2000Vector(lookB_copy); // Get body look vector in J2000 (Earth) frame
  std::vector<double> lookC = m_camera->instrumentRotation()->ReferenceVector(lookJ); // Get J2000 look vector in camera frame
  
  // Convert vector format
  Vector3 lookI;
  std::copy( lookC.begin(), lookC.end(), lookI.begin() );

  //std::cout << "lookI = " << lookI << std::endl;
  
  // Apply extra rotation to the look vector (equal to the one in the function below)
  // R_offset = instrument(good)_from_instrument
  vw::math::Matrix<double,3,3> R_offset = vw::math::euler_to_rotation_matrix(m_rotation[0], m_rotation[1], m_rotation[2], "xyz");
  Vector3 look = R_offset * lookI;
  
  //std::cout << "look = " << look << std::endl;
  
  // Projecting to mm focal plane
  look = m_camera->FocalLength() * (look / look[2]);
  m_distortmap->SetUndistortedFocalPlane(look[0], look[1]);
  m_focalmap->SetFocalPlane( m_distortmap->FocalPlaneX(),
                             m_distortmap->FocalPlaneY() );
  result_type result(1);
  // Not exactly sure about lineoffset .. but ISIS does it
  result[0] = m_focalmap->DetectorLineOffset() - m_focalmap->DetectorLine();

  return result;
}

Vector2
IsisInterfaceLineScanRot::point_to_pixel_rotated( vw::Vector3 const& point, vw::Vector3 const& rotAngles) const 
{

  // First seed LMA with an ephemeris time in the middle of the image
  double middle = lines() / 2;
  m_detectmap->SetParent( 1, m_alphacube.AlphaLine(middle) );
  double start_e = m_camera->time().Et();

  // Build LMA
  EphemerisLMA_rot model( point, rotAngles, m_camera.get(), m_distortmap, m_focalmap );
  int status;
  Vector<double> objective(1), start(1);
  start[0] = start_e;
  Vector<double> solution_e = math::levenberg_marquardt( model,
                                                         start,
                                                         objective,
                                                         status );

  // Make sure we found ideal time
  VW_ASSERT( status > 0,
             MathErr() << " Unable to project point into linescan camera " );

  // Converting now to pixel
  m_camera->setTime(Isis::iTime( solution_e[0] ));

  // Working out pointing
  m_camera->instrumentPosition(&m_center[0]); // Get camera position in GCC coordinates
  m_center *= 1000; // Convert to meters
  Vector3 look = normalize(point-m_center); // Vector from point to camera in GCC coordinates

  // Calculating Rotation to camera frame
  std::vector<double> rot_inst = m_camera->instrumentRotation()->Matrix();
  std::vector<double> rot_body = m_camera->bodyRotation()->Matrix();
  MatrixProxy<double,3,3> R_inst(&(rot_inst[0])); // Rotation of spacecraft and instrument relative to J2000
  MatrixProxy<double,3,3> R_body(&(rot_body[0])); // Rotation of planet relative to J2000 (earth-centered) frame
  
  //DEBUG: Extract instrument rotation angles
  Vector3 instAngles = vw::math::rotation_matrix_to_euler_xyz(R_inst);
  
  //std::cout << "InstAngles = " << instAngles << std::endl;
  double rad2deg = 180/M_PI;
  //printf("InstAngles (deg) = %lf, %lf, %lf\n", instAngles[0]*rad2deg, instAngles[1]*rad2deg, instAngles[2]*rad2deg);
  
  Vector3 bAngles = vw::math::rotation_matrix_to_euler_xyz(R_body);
  //printf("Body Angles (deg) = %lf, %lf, %lf\n", bAngles[0]*rad2deg, bAngles[1]*rad2deg, bAngles[2]*rad2deg);
  
  
  // Add in input offsets
  Vector3 offsetAngles = /*instAngles +*/ rotAngles;
  vw::math::Matrix<double,3,3> R_offset = vw::math::euler_to_rotation_matrix(offsetAngles[0], offsetAngles[1], offsetAngles[2], "xyz");
  // Try partial rotation solution.
  //vw::math::Matrix<double,3,3> R_offset = vw::math::euler_to_rotation_matrix(offsetAngles[0], offsetAngles[1], 0, "xyz");
    
  // Existing offset matrix (from the fk file) = instrument(bad)_from_sc
  // New offset matrix (from input params) = instrument(good)_from_instrument(bad)
  
  // Generate rotation matrix from input angles
  //vw::math::Quat offsetRot(rotAngles);
  //vw::math::Matrix<double,3,3> R_offset;
  //offsetRot.rotation_matrix(R_offset);
  
  // Compute pose --> rotation from camera coordinates to body coordinates ?
  //m_pose = Quat(R_body*transpose(R_inst)); // Original
  //m_pose = Quat(R_body * transpose(R_offset * R_inst)); // Trial 1
  //m_pose = Quat(R_body * transpose(R_inst * R_offset)); // Trial 3 = fail
  //m_pose = Quat(R_body * transpose(R_inst) * R_offset); // Trial 4
  //m_pose = Quat(R_body * R_offset * transpose(R_inst)); // Trial 6
  //m_pose = Quat(R_body * transpose(R_inst * R_offset)); // Trial 8
  //m_pose = Quat(R_offset * R_body*transpose(R_inst)); // Trial 10
  
  
  m_pose = Quat(R_body*transpose(R_offset * R_inst)); // Trial 12
  
  //std::cout << "R_body   = " << R_body << std::endl;
  //std::cout << "R_inst   = " << R_inst << std::endl;
  //std::cout << "R_offset = " << R_offset << std::endl;
  //std::cout << "R_offset*R_inst = " << R_offset*R_inst << std::endl;
  //std::cout << "pose     = " << m_pose << std::endl;
  //vw::math::Matrix<double,3,3> R_final = R_body*transpose(R_offset * R_inst);
  
  // R_body = body_from_J2000
  // R_inst = instrument_from_J2000 = instrument_from_sc * sc_from_J2000
  // transpose(R_inst) = J2000_from_instrument
  
  // R_body*transpose(R_inst) = body_from_J2000 * J2000_from_instrument = body_from_instrument
  // R_offset = instrument(good)_from_instrument
  // transpose(R_offset * R_inst) = J2000_from_instrument(good)
  
  // m_pose = body_from_instrument
  // inverse(m_pose) = instrument_from_body
  
  look = inverse(m_pose).rotate( look ); // Apply rotation from body coordinates to camera coordinates to look vector
  look = m_camera->FocalLength() * ( look / look[2] );
  m_distortmap->SetUndistortedFocalPlane( look[0], look[1] );
  m_focalmap->SetFocalPlane( m_distortmap->FocalPlaneX(),
                             m_distortmap->FocalPlaneY() );
  m_detectmap->SetDetector( m_focalmap->DetectorSample(),
                            m_focalmap->DetectorLine() );
  Vector2 pixel( m_detectmap->ParentSample(),
                 m_detectmap->ParentLine() );
  pixel[0] = m_alphacube.BetaSample( pixel[0] );
  pixel[1] = m_alphacube.BetaLine( pixel[1] );
  SetTime( pixel, false );

  pixel -= Vector2(1,1);
  return pixel;
}

} // end namespace isis
} // end namespace asp



//----------------------------------------------------------------------------------------------------------------



struct Parameters : asp::BaseOptions 
{
  // Input paths
  std::string leftFilePath;
  std::string rightFilePath;
  std::string outputPath;
  std::string gdcPointsOutPath;
  std::string matchingPointsPath; ///< If not set, ipfind/ransac is used to find matching points

  bool worldTransform;  ///< Solve for transform in world coordinate frame instead of camera frame
  bool includePosition; ///< Also solve for a position offset in the specified coordinate frame

  int  cropWidth; ///< Specifies image overlap for use with ipfind
};



bool handle_arguments(int argc, char* argv[],
                     Parameters &opt) 
{ 
  po::options_description general_options("Options");
  general_options.add_options()
    ("outputPath",         po::value      (&opt.outputPath        )->default_value(""),     "Write angles to this path (in degrees)")
    ("matchingPixelsPath", po::value      (&opt.matchingPointsPath)->default_value(""),     "File to load matching points from")
    ("gdcPointsOutPath",   po::value      (&opt.gdcPointsOutPath  )->default_value(""),     "Output path for final GDC points")
    ("worldTransform",     po::bool_switch(&opt.worldTransform    )->default_value(false),  "Compute transform in world frame instead of camera frame")
    ("includePosition",    po::bool_switch(&opt.includePosition   )->default_value(false),  "Also solve for a 3D position offset (in chosen frame)")
    ("crop-width",         po::value      (&opt.cropWidth         )->default_value(200),    "Crop images to this width before disparity search");
  
  general_options.add( asp::BaseOptionsDescription(opt) );
    
  po::options_description positional("");
  positional.add_options()
    ("left",  po::value(&opt.leftFilePath))
    ("right", po::value(&opt.rightFilePath));  
    
  po::positional_options_description positional_desc;
  positional_desc.add("left",  1);
  positional_desc.add("right", 1);
  

  std::string usage("[options] <left> <right>");
  po::variables_map vm =
    asp::check_command_line( argc, argv, opt, general_options, general_options,
                             positional, positional_desc, usage );

  if ( !vm.count("left") || !vm.count("right") )
    vw_throw( ArgumentErr() << "Requires <left> and <right> input in order to proceed.\n\n"
              << usage << general_options );

  return true;
}

//----------------------------------------------------------------------------------------------------------------
/// Class for solving for the rotation between two LRONAC cameras
class LrocPairModel : public vw::math::LeastSquaresModelBase<LrocPairModel>
{
public:  // Definitions
  
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
  
private: // Variables
  
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
  
public: // Functions  
  
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
                               Vector<double> &stateEstimate, Vector<double> &packedObsVector)
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

    //DEBUG
    // Set up georeference class with default moon datum
    vw::cartography::Datum datum("D_MOON");
  
  
    printf("Setting up state estimate\n");

    // For each input point pair
    Vector3 pointLoc, lastPointLoc;
    for (size_t i=0; i<numPoints; ++i)
    {
      // Set up point pair
      Vector2 leftPixel (leftCols [i], leftRows [i]);
      Vector2 rightPixel(rightCols[i], rightRows[i]);
      
      // Compute the intersection location
      Vector3 zeroRotVec(0, 0, 0);
      //_rightCameraRotatedModel.set_axis_angle_rotation(zeroRotVec);
      double triangulationError;
      pointLoc = _stereoModel(leftPixel, rightPixel, triangulationError); // Not working!
      
      Vector3 leftCamCenter  = _leftCameraModel.camera_center(leftPixel);
      Vector3 rightCamCenter = _rightCameraModel.camera_center(rightPixel);
      
      Quaternion<double> leftCamPose  = vw::math::normalize(_leftCameraModel.camera_pose (leftPixel));
      Quaternion<double> rightCamPose = vw::math::normalize(_rightCameraModel.camera_pose(rightPixel));
      Vector3            leftVec      = vw::math::normalize(_leftCameraModel.pixel_to_vector (leftPixel));
      Vector3            rightVec     = vw::math::normalize(_rightCameraModel.pixel_to_vector(rightPixel));
      //double             vectorAngle  = acos(vw::math::dot_prod(leftVec, rightVec));
      Quaternion<double> conjLeftPose = vw::math::conj(leftCamPose);
      Quaternion<double> invLeftPose  = vw::math::inverse(leftCamPose);
      Quaternion<double> rotDiff      = rightCamPose*invLeftPose;
      Vector3 axis; double angle;
      rotDiff.axis_angle(axis, angle);

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
  Vector3 midPoint = 0.5 * (closestPoint1 + closestPoint2);
  //pointLoc = midPoint; // HIJACK TRIANGULATION CALCULATIONS!
      
  Vector3 errorVec = closestPoint1 - closestPoint2;
  triangulationError = vw::math::norm_2(errorVec);
          
      //printf("\n");
      //std::cout << "closestPoint1   = " << closestPoint1 <<", radius = " << vw::math::norm_2(closestPoint1)/1000.0 << std::endl;
      //std::cout << "closestPoint2   = " << closestPoint2 <<", radius = " << vw::math::norm_2(closestPoint2)/1000.0 << std::endl;
      //std::cout << "errorVec        = " << errorVec      << std::endl;
      //std::cout << "projection dist = " << vw::math::norm_2(closestPoint1 - leftCamCenter)/1000.0 << std::endl;
          
      
      //std::cout.precision(16);
      //std::cout << "Left  camera center = " << leftCamCenter [0] <<", "<< leftCamCenter [1] <<", "<< leftCamCenter [2]<<", radius = " << vw::math::norm_2(leftCamCenter)/1000.0 << std::endl;
      //std::cout << "Right camera center = " << rightCamCenter[0] <<", "<< rightCamCenter[1] <<", "<< rightCamCenter[2]<<", radius = " << vw::math::norm_2(rightCamCenter)/1000.0 << std::endl;
      //std::cout << "Center diff         = " << leftCamCenter[0]-rightCamCenter[0] <<", "<< leftCamCenter[1]-rightCamCenter[1] <<", "<< leftCamCenter[2]-rightCamCenter[2] << std::endl;
      //printf("Angle between cameras = %lf\n", angle);
      //printf("Center abs diff = %lf\n", sqrt( pow(leftCamCenter[0]-rightCamCenter[0], 2) + pow(leftCamCenter[1]-rightCamCenter[1], 2) + pow(leftCamCenter[2]-rightCamCenter[2], 2) ));
      
      //printf("Desired abs diff = %lf\n", sqrt( pow(134.62-101.60, 2) + pow(88.90-88.90, 2) + pow(-17.78 - -17.78, 2) ));
      // = 0.33 meters, from the kernel definition file
      // Angular difference on start should be over 2.5 degrees
      
      std::cout << "Left pixel  = " << leftPixel << " Right pixel = " << rightPixel << std::endl;
      std::cout << "Initial point " << i << " = " << pointLoc << " Triangulation error = " << triangulationError << std::endl;
      
      
      Vector2 rightProjection, leftProjection;
      /*
      if (_solveWorldFrame)
        rightProjection = _rightCameraRotatedModel.point_to_pixel(pointLoc);
      else // Camera frame
        //rightProjection = _rightCameraModel.point_to_pixel_rotated(pointLoc, Vector3(0,0,0));
        rightProjection = _rightCameraModel.point_to_pixel_rotated(pointLoc, Vector3(0,0,0));
      */
      leftProjection  = _leftCameraModel.point_to_pixel(pointLoc);
      rightProjection = _rightCameraModel.point_to_pixel(pointLoc);
      std::cout << "**Left: "  << leftProjection  << " - " << leftPixel  << std::endl;
      std::cout << "**Right: " << rightProjection << " - " << rightPixel << std::endl;
      
      // Sanity check
      const double intersectionRadius = vw::math::norm_2(pointLoc);
      const double camRadius          = vw::math::norm_2(rightCamCenter);
      //std::cout << "intersection radius = " << intersectionRadius << std::endl;
      
      // Convert from GCC to GDC
      Vector3 gdcCoord = datum.cartesian_to_geodetic(midPoint);
      
      std::cout << "elevation above ellipsoid = " << gdcCoord[2] << std::endl;
      
      //// Try dropping the elevation down to the datum level
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
        rightProjection = _rightCameraRotatedModel.point_to_pixel(thisPoint);
      else // Camera frame
        rightProjection = _rightCameraModel.point_to_pixel_rotated(thisPoint, rotVec);
      
                        
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
  
  
  /// * The domain_type must implement a method: domain_type domain_type::operator+( gradient_type const& g ) const;
  ///   that adds a tangent vector to a position.  You get this for
  ///   free if both domain_type and gradient_type are Vector<foo>.
  ///   This is where you do most of the hard work if domain_type
  ///   represents some non-trivial topological space.
  
  
  /// * The result_type must implement a method: double result_type::norm_2( result_type const& g ) const;
  ///   that is used in some optimizers to compute the error.  You get this for
  ///   free if result_type is a Vector<foo>.
  
  
  /// * The jacobian_type must implement several matrix-like methods such as
  ///   scalar multiplication on the left.  You get this for free in the usual case when
  ///   jacobian_type is just Matrix<foo>.
  ///
  
}; // End class LrocPairModel
//----------------------------------------------------------------------------------------------------------------

// Load mathing points from a file
bool loadMatchingPixels(const std::string &pointPath,
          Vector<double> &leftRow, Vector<double> &leftCol, Vector<double> &rightRow, Vector<double> &rightCol)
{
  // The input file is a line for each point: sample1, line1, sample2, line2
  std::ifstream file(pointPath.c_str());
  if (file.fail())
  {
    printf("Failed to open point file %s\n", pointPath.c_str());
    return false;
  }
  
  char comma;
  std::string line;
  double sample1, line1, sample2, line2;
  std::vector<double> row1, row2, col1, col2;
  row1.reserve(50);
  row2.reserve(50);
  col1.reserve(50);
  col2.reserve(50);
  int p=0;
  while (std::getline(file, line)) // Read in each line and append to vectors
  {
    if (line.size() < 8) // Stop if we hit a blank line
      break;
    std::stringstream s(line);
    s >> sample1 >> comma >> line1 >> comma >> sample2 >> comma >> line2;
    col1.push_back(sample1);
    row1.push_back(line1);
    col2.push_back(sample2);
    row2.push_back(line2);
    printf("p: %d, --> %lf, %lf, %lf, %lf\n", p, sample1, line1, sample2, line2);
    ++p;
  }
  file.close();
  
  leftRow.set_size (p);
  leftCol.set_size (p);
  rightRow.set_size(p);
  rightCol.set_size(p);
  for (int i=0; i<p; ++i)
  {
    leftCol [i] = col1[i];
    leftRow [i] = row1[i];
    rightCol[i] = col2[i];
    rightRow[i] = row2[i];
  }
  
  return true;
}

// Search for matching pixels in the LE/RE overlap
bool findMatchingPixels(const Parameters &params, Vector<double> &leftRow, Vector<double> &leftCol, Vector<double> &rightRow, Vector<double> &rightCol)
{
 
  // Load both images  
  printf("Loading images left=%s and right=%s...\n",
         params.leftFilePath.c_str(),
         params.rightFilePath.c_str());
  DiskImageView<PixelGray<float> > left_disk_image (params.leftFilePath );
  DiskImageView<PixelGray<float> > right_disk_image(params.rightFilePath);
  
  printf("Left  input image size: %d rows, %d cols\n", left_disk_image.rows(),  left_disk_image.cols());
  printf("Right input image size: %d rows, %d cols\n", right_disk_image.rows(), right_disk_image.cols());
  
  //const int imageWidth      = std::min(left_disk_image.cols(), right_disk_image.cols());
  const int imageHeight     = std::min(left_disk_image.rows(), right_disk_image.rows());
  const int imageTopRow     = 0;
  //const int imageMidPointX  = imageWidth / 2;
  //const int cropStartX      = imageMidPointX - (params.cropWidth/2);

  // Restrict processing to the border of the images
  // - Since both images were nproj'd the overlap areas should be in about the same spots.
  const int leftStartX = left_disk_image.cols()-params.cropWidth;
  const BBox2i leftRoi (leftStartX, imageTopRow, params.cropWidth, imageHeight);
  const BBox2i rightRoi(0,          imageTopRow, params.cropWidth, imageHeight);
  std::cout << "Left  overlap ROI = " << leftRoi  << std::endl;
  std::cout << "Right overlap ROI = " << rightRoi << std::endl;


  // Now use interest point finding/matching functions to estimate the search offset between the images
  printf("Gathering interest points...\n");

  // Gather interest points
  asp::IntegralAutoGainDetector detector( 500 );
  ip::InterestPointList ip1 = ip::detect_interest_points( vw::create_mask_less_or_equal(crop(left_disk_image,  leftRoi ), 0), detector );
  ip::InterestPointList ip2 = ip::detect_interest_points( vw::create_mask_less_or_equal(crop(right_disk_image, rightRoi), 0), detector );
  printf("Found %lu, %lu interest points.\n", ip1.size(), ip2.size());
      
  ip::SGradDescriptorGenerator descriptor;
  describe_interest_points( vw::create_mask_less_or_equal(crop(left_disk_image,  leftRoi ), 0), descriptor, ip1 );
  describe_interest_points( vw::create_mask_less_or_equal(crop(right_disk_image, rightRoi), 0), descriptor, ip2 );

  // Match interest points
  ip::DefaultMatcher matcher(0.5);
  std::vector<ip::InterestPoint> matched_ip1, matched_ip2;
  matcher(ip1, ip2, matched_ip1, matched_ip2 );
  ip::remove_duplicates( matched_ip1, matched_ip2 );

  if (matched_ip1.empty() || matched_ip2.empty())
  {
    printf("Failed to find any matching interest points, defaulting to large search range.\n");
    return false;
  }
  
  printf("Found %lu matched interest points.\n", matched_ip1.size());

  // Init inlier indices to all matched interest points
  std::vector<size_t> inlierIndices(matched_ip1.size());
  for (size_t i=0; i<matched_ip1.size(); ++i)
    inlierIndices[i] = i;

  printf("Filtering points with RANSAC...\n");
  // Filter interest point matches
  int    numIterations       = 100;
  double inlierThreshold     = 5.0; // Want to be somewhat generous here
  int    minNumOutputInliers = 100;
  math::RandomSampleConsensus<math::SimilarityFittingFunctor, math::InterestPointErrorMetric> ransac( math::SimilarityFittingFunctor(),
                          math::InterestPointErrorMetric(),
                          numIterations, inlierThreshold, minNumOutputInliers, false );
  std::vector<Vector3> ransac_ip1 = ip::iplist_to_vectorlist(matched_ip1);
  std::vector<Vector3> ransac_ip2 = ip::iplist_to_vectorlist(matched_ip2);

  // Try to find a consistent offset between the matched interest points
  try
  {
    // Find best transform
    Matrix<double> H(ransac(ransac_ip1, ransac_ip2));
    std::cout << "ipfind based similarity: " << H << std::endl;

    // Get list of interest points consistent with the transform
    std::vector<size_t> inlierIndicesF = ransac.inlier_indices(H, ransac_ip1, ransac_ip2); //TODO: Why is this failing?
    size_t numInliers = inlierIndicesF.size(); 
    printf("Found %lu inliers\n", numInliers);
    
    std::vector<Vector3> tempL(numInliers), tempR(numInliers);
    for (size_t i=0; i<numInliers; ++i)  // Replace contents of ransac_ip1 with inliers only
    {
      size_t index = inlierIndicesF[i];
      tempL[i] = ransac_ip1[index];
      tempR[i] = ransac_ip2[index];
    }
    
    ransac_ip1 = tempL;
    ransac_ip2 = tempR;        
  }
  catch(...) // Handle a RANSAC failure
  {
    printf("RANSAC solution failed!\n");
    return false;
  }

  // Now that we have correspondence points, feed them into an angular solver.

  // Convert the matching points into the correct format!
  const int    pointSkip     = 500; //TODO: MOVE THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  const size_t numMatchedPts = ransac_ip1.size() / pointSkip;
  printf("Num sampled points = %lu\n", numMatchedPts);
  leftRow.set_size(numMatchedPts), leftCol.set_size(numMatchedPts), rightRow.set_size(numMatchedPts), rightCol.set_size(numMatchedPts);
  int i = 0;
  for (size_t p=0; p<numMatchedPts; ++p)
  {
    leftCol [p] = ransac_ip1[i][0] + leftStartX;
    leftRow [p] = ransac_ip1[i][1];
    rightCol[p] = ransac_ip2[i][0];
    rightRow[p] = ransac_ip2[i][1];
    //printf("p: %d, i: %d --> %lf, %lf, %lf, %lf\n", p, i, leftCol[p], leftRow[p], rightCol[p], rightRow[p]);
    i+= pointSkip;
  }  
  
  return true;
}


// Main solver function
bool optimizeRotations(Parameters & params)
{
  
  // Verify images are present
  boost::filesystem::path leftBoostPath (params.leftFilePath );
  boost::filesystem::path rightBoostPath(params.rightFilePath);
  
  if (!boost::filesystem::exists(boost::filesystem::path(params.leftFilePath)))
  {
    printf("Error: input file %s is missing!\n", params.leftFilePath.c_str());
    return false;
  }
  if (!boost::filesystem::exists(boost::filesystem::path(params.rightFilePath)))
  {
    printf("Error: input file %s is missing!\n", params.rightFilePath.c_str());
    return false;
  }

  printf("Constructing geometry class\n");

  // Initialize the geometry/solver class for the two input cubes
  LrocPairModel lrocClass(params.leftFilePath, params.rightFilePath, params.worldTransform, params.includePosition);
  
  //DEBUG - call functions with fixed points to print some pixel info (no rotation)
  lrocClass.getRightPixelRot(Vector3(-1.09977e+06, 387050, 1.29165e+06), Vector3());
  lrocClass.getRightPixelRot(Vector3(-1.09272e+06, 384682, 1.29613e+06), Vector3());
  lrocClass.getRightPixelRot(Vector3(-1.08521e+06, 382190, 1.2999e+06 ), Vector3());
  lrocClass.getRightPixelRot(Vector3(-1.07943e+06, 380243, 1.30479e+06), Vector3());

  
  Vector<double> leftRow, leftCol, rightRow, rightCol;
  
  if (params.matchingPointsPath.empty())
  {
    printf("Searching for matching pixels in image overlap region\n");
    if (!findMatchingPixels(params, leftRow, leftCol, rightRow, rightCol))
      return false;
  }
  else
  {
    printf("Loading list of matched pixels from file\n");
    if (!loadMatchingPixels(params.matchingPointsPath, leftRow, leftCol, rightRow, rightCol))
      return false;
  }


  const size_t numMatchedPts = leftRow.size();

  // Load the inital points into the solver
  printf("Initializing solver state...\n");
  Vector<double> initialState, packedObservations;
  if (!lrocClass.getInitialStateEstimate(leftRow, leftCol, rightRow, rightCol, initialState, packedObservations))
  {
    printf("Failed to get initial state!\n");
    return false;
  }


  //TODO: Move these
  std::string initialErrorPath    = "/home/smcmich1/initialAngleError.csv";
  std::string finalErrorPath      = "/home/smcmich1/finalAngleError.csv";
  std::string initialStatePath    = "/home/smcmich1/initialAngleState.csv";
  std::string finalStatePath      = "/home/smcmich1/finalAngleState.csv";
  std::string finalStateDiffPath  = "/home/smcmich1/finalAngleStateDiff.csv";
  std::string initialGdcCoordPath = "/home/smcmich1/initialGdcPoints.csv";
  std::string finalGdcCoordPath   = "/home/smcmich1/finalGdcPoints.csv";
  
  // Set up georeference class with default moon datum
  vw::cartography::Datum datum("D_MOON");
  
  //printf("Writing initial state log to %s\n", initialStatePath.c_str());
  std::ofstream initialObsFile("/home/smcmich1/initialObsVector.csv"); 
  for (size_t i=0; i<packedObservations.size(); ++i)
    initialObsFile << packedObservations[i] << std::endl;
  initialObsFile.close();
  
  //return false;
  
  Vector<double> initialPredictions = lrocClass(initialState);
  std::ofstream initialPredFile("/home/smcmich1/initialObsPrediction.csv"); 
  for (size_t i=0; i<initialPredictions.size(); ++i)
    initialPredFile << initialPredictions[i] << std::endl;
  initialPredFile.close();
   
  
  //printf("Writing initial state log to %s\n", initialStatePath.c_str());
  std::ofstream initialStateFile(initialStatePath.c_str()); 
  for (size_t i=0; i<initialState.size(); ++i)
    initialStateFile << initialState[i] << std::endl;
  initialStateFile.close();
  
  int startOfPts = 3;
  if (params.includePosition)
    startOfPts = 6;
  
  // Write initial points as GDC coordinates for google earth
  std::ofstream initialGdcCoordFile(initialGdcCoordPath.c_str());
  for (size_t i=0; i<numMatchedPts; ++i)
  {
    // Convert from GCC to GDC
    Vector3 gccPoint(initialState[startOfPts+ 3*i], initialState[startOfPts+ 3*i+1], initialState[startOfPts+ 3*i+2]);
    Vector3 gdcCoord = datum.cartesian_to_geodetic(gccPoint);
    initialGdcCoordFile.precision(12);
    initialGdcCoordFile << gdcCoord[0] << std::endl;
    initialGdcCoordFile << gdcCoord[1] << std::endl;
    initialGdcCoordFile << gdcCoord[2] << std::endl;
  }
  initialGdcCoordFile.close();
  
  // Compute the initial error - euclidean point distance
  //printf("Writing initial error log to %s\n", initialErrorPath.c_str());
  std::ofstream initialErrorFile(initialErrorPath.c_str()); 
  double meanInitialError = 0;  
  std::vector<double> currentError = lrocClass.computeError(initialState);
  for (size_t i=0; i<currentError.size(); ++i)
  {   
    initialErrorFile << currentError[i] << std::endl;
    meanInitialError += currentError[i];
  }
  initialErrorFile.close();
  meanInitialError = meanInitialError / currentError.size();
  printf("Mean point error before optimization = %lf\n", meanInitialError);
  
  Vector<double> initialComputedObservations = lrocClass(initialState);

// DEBUG: Check output Jacobian
  Matrix<double> initialJac = lrocClass.jacobian(initialState);
  std::ofstream initialJacFile("/home/smcmich1/initialJac.csv");
  initialJacFile.precision(3);
  for (size_t r=0; r<initialJac.rows(); ++r)
  {
    for (size_t c=0; c<initialJac.cols(); ++c)
    {
      initialJacFile << initialJac(r,c) << " ";
    }
    initialJacFile << std::endl;
  }
  initialJacFile.close();
  
  
  printf("Running solver...\n");

  // Now pass geometry class into the solver function
  int status;
  Vector<double> finalParams = vw::math::levenberg_marquardt(lrocClass, initialState, packedObservations, status);
//                                                   double abs_tolerance = VW_MATH_LM_ABS_TOL,
//                                                   double rel_tolerance = VW_MATH_LM_REL_TOL,
//                                                   double max_iterations = VW_MATH_LM_MAX_ITER)
  std::cout << "Status = " << status << std::endl;

  
  //printf("Writing final error log to %s\n", finalErrorPath.c_str());
  std::ofstream finalErrorFile(finalErrorPath.c_str()); 
  double meanFinalError = 0;
  std::vector<double> finalError = lrocClass.computeError(finalParams);
  Vector<double> finalPredictions = lrocClass(finalParams);
  Vector<double> rawError(finalPredictions.size());
  std::ofstream predictionFile("/home/smcmich1/finalPredictions.csv");
  for (size_t i=0; i<finalPredictions.size(); ++i)
  {
    predictionFile << finalPredictions[i] << std::endl;
    rawError[i] = packedObservations[i] - finalPredictions[i];
    //finalErrorFile << packedObservations[i] - finalPredictions[i] << std::endl;
    if (i % 4 == 0)
      meanFinalError += finalError[i/4];
  }
  
  for (size_t i=0; i<finalError.size(); ++i)
  {   
    finalErrorFile << finalError[i] << std::endl;
    meanFinalError += finalError[i];
  }
  meanFinalError = meanFinalError / finalError.size(); 
  
  
  predictionFile.close();
  finalErrorFile.close();
  //meanFinalError = meanFinalError / currentError.size(); 
  
  
  
  
  
  // DEBUG: Check output Jacobian -------------------------------------------
  Matrix<double> finalJac = lrocClass.jacobian(finalParams);
  std::ofstream finalJacFile("/home/smcmich1/finalJac.csv");
  finalJacFile.precision(3);
  for (size_t r=0; r<finalJac.rows(); ++r)
  {
    for (size_t c=0; c<finalJac.cols(); ++c)
    {
      finalJacFile << finalJac(r,c) << " ";
    }
    finalJacFile << std::endl;
  }
  finalJacFile.close();
  
  
  
  //printf("Writing final state log to %s\n", finalStatePath.c_str());
  std::ofstream finalStateFile(finalStatePath.c_str()); 
  for (size_t i=0; i<finalParams.size(); ++i)
    finalStateFile << finalParams[i] << std::endl;
  finalStateFile.close();

  // Write initial points as GDC coordinates for google earth
  std::ofstream finalGdcCoordFile;
  if (params.gdcPointsOutPath.empty())
    finalGdcCoordFile.open(finalGdcCoordPath.c_str()); // TODO: Remove deubg default
  else
    finalGdcCoordFile.open(params.gdcPointsOutPath.c_str());
  for (size_t i=0; i<numMatchedPts; ++i)
  {
    // Convert from GCC to GDC
    Vector3 gccPoint(finalParams[startOfPts+ 3*i], finalParams[startOfPts+ 3*i+1], finalParams[startOfPts+ 3*i+2]);
    Vector3 gdcCoord = datum.cartesian_to_geodetic(gccPoint);
    finalGdcCoordFile.precision(12);
    finalGdcCoordFile << gdcCoord[0] << std::endl;
    finalGdcCoordFile << gdcCoord[1] << std::endl;
    finalGdcCoordFile << gdcCoord[2] << std::endl;
  }
  finalGdcCoordFile.close();

  std::ofstream finalStateDiffFile(finalStateDiffPath.c_str()); 
  for (size_t i=0; i<finalParams.size(); ++i)
    finalStateDiffFile << finalParams[i] - initialState[i] << std::endl;
  finalStateDiffFile.close();
  
  printf("Mean point error after optimization = %lf\n", meanFinalError);
  printf("Mean error change = %lf\n", meanFinalError - meanInitialError);
  
  const double rad2deg = 180.0 / M_PI;
  const double deg2rad = M_PI / 180.0;
  printf("Adjustment rotation angles (degrees): %lf, %lf, %lf\n", finalParams[0]*rad2deg, finalParams[1]*rad2deg, finalParams[2]*rad2deg);
  //printf("Output rotation angles (degrees): %lf, %lf, %lf\n", finalParams[0]*rad2deg, finalParams[1]*rad2deg, finalParams[2]*rad2deg);
  if (params.includePosition)
    printf("Adjustment offset (meters): %lf, %lf, %lf\n", finalParams[3], finalParams[4], finalParams[5]);
  
  
  // These are instrument rotations from the frame file - rotation order is R = Rx(1.29)*Ry(-0.079)*Rz(180)
  // sc_from_instrument,  ( 1.129, -0.079, 180.0 ) [units in degrees]
  
  // vw::math::euler_to_rotation_matrix(a, b, c, 'order') --> R(c)*R(b)*R(a), order specifies which rotation axes
  
  //vw::math::Matrix<double,3,3> frameRotation = vw::math::euler_to_rotation_matrix(180.0*deg2rad, -0.079*deg2rad, 1.129*deg2rad, "zyx");
  vw::math::Matrix<double,3,3> frameRotation = vw::math::euler_to_rotation_matrix(1.129*deg2rad, -0.079*deg2rad, 180.0*deg2rad, "xyz");
   std:: cout << "Base rotation angles (degrees) = Rx(1.129) * Ry(-0.079) * Rz(180.0)"  << std::endl;

  // Existing offset matrix (from the fk file) = sc_from_instrument
  
  // This is the rotation matrix we solved for (matches code in IsisInterfaceLineScan.cc)
  // solvedRotation = instrument(good)_from_instrument
  vw::math::Matrix<double,3,3> solvedRotation = vw::math::euler_to_rotation_matrix(finalParams[0], finalParams[1], finalParams[2], "xyz");
  std::cout << "Solved rotation: " << solvedRotation << std::endl;
  // Multiply them to get the net rotation

  //vw::math::Matrix<double,3,3> outputRotation = frameRotation; //TEST: Make sure original rotations kept intact!
  //vw::math::Matrix<double,3,3> fixed_from_sc  = transpose(solvedRotation) * frameRotation;
  //vw::math::Matrix<double,3,3> outputRotation = transpose(fixed_from_sc); // == sc_from_fixed
  
  vw::math::Matrix<double,3,3> outputRotation;
  if (params.worldTransform)
    outputRotation = solvedRotation;
  else
    outputRotation = transpose(frameRotation * transpose(solvedRotation)); // == sc_from_instrument(good)
  std::cout << "Output rotation: " << outputRotation << std::endl;
  
  //vw::math::Matrix<double,3,3> outputRotationT = transpose(frameRotation * transpose(solvedRotation));
  //std::cout << "Output rotation TRANS: " << outputRotationT << std::endl;
  
  // Extract rotation angles --> Returns Rx(a)*Ry(b)*Rz(c)
  // -- Angle order needs to be  1, 2, 3
  //Vector3 outputAngles = vw::math::rotation_matrix_to_euler_zyx(outputRotation); //bad
  Vector3 outputAngles = vw::math::rotation_matrix_to_euler_xyz(outputRotation); // good?
  // Output angles (rZ, rY, rX)
  
  printf("Combined rotation angles (degrees): %lf, %lf, %lf\n", outputAngles[2]*rad2deg, outputAngles[1]*rad2deg, outputAngles[0]*rad2deg);
  
  if (!params.outputPath.empty()) // Dump rotation angles to a simple text file (in degrees)
  {
    printf("Writing output file %s\n", params.outputPath.c_str()); 
    std::ofstream outputFile(params.outputPath.c_str());
    //outputFile << outputAngles[0]*rad2deg << std::endl; // X rotation
    //outputFile << outputAngles[1]*rad2deg << std::endl; // Y rotation
    //outputFile << outputAngles[2]*rad2deg << std::endl; // Z rotation

    // Write out rotation matrix
    for (unsigned int r=0; r<outputRotation.rows(); ++r)
    {
      for (unsigned int c=0; c<outputRotation.cols(); ++c)
      {
        outputFile << outputRotation(r,c) << std::endl;
      }
    }
    
    if (params.includePosition) // Also write out the translation
      outputFile << finalParams[3] << finalParams[4] << finalParams[5] << std::endl;
    
    outputFile.close();
  }

  //DEBUG - call functions with fixed points to print some pixel info (solved rotation)
  lrocClass.getRightPixelRot(Vector3(-1.09977e+06, 387050, 1.29165e+06), Vector3(finalParams[0], finalParams[1], finalParams[2]));
  lrocClass.getRightPixelRot(Vector3(-1.09272e+06, 384682, 1.29613e+06), Vector3(finalParams[0], finalParams[1], finalParams[2]));
  lrocClass.getRightPixelRot(Vector3(-1.08521e+06, 382190, 1.2999e+06 ), Vector3(finalParams[0], finalParams[1], finalParams[2]));
  lrocClass.getRightPixelRot(Vector3(-1.07943e+06, 380243, 1.30479e+06), Vector3(finalParams[0], finalParams[1], finalParams[2]));

  return true;
}

int main(int argc, char* argv[]) 
{ 
  try 
  {
    // Parse the input parameters
    Parameters params;
    if (!handle_arguments(argc, argv, params))
    {
      printf("Failed to parse input parameters!\n");
      return false;
    }
    
    optimizeRotations(params);
  } ASP_STANDARD_CATCHES;

  return 0;
}






