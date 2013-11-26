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

#ifndef LRONAC_SOLVER_SUPPORT_H
#define LRONAC_SOLVER_SUPPORT_H


/// \file lronacAngleSolverSupport.cc
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

using namespace vw;
using namespace vw::stereo;
using namespace asp::isis;
using namespace asp;
using std::endl;
using std::setprecision;
using std::setw;



//TODO: Replace this code!!!!


bool raySphereIntersect(const Vector3 &rayOrigin, const Vector3 &rayVector, const double radius, Vector3 &intersection)
{
    // Get two points on the ray
    double xA = rayOrigin[0]; 
    double yA = rayOrigin[1];
    double zA = rayOrigin[2];

    double a = rayVector[0]*rayVector[0] + rayVector[1]*rayVector[1] + rayVector[2]*rayVector[2];
    double b = 2*(rayVector[0]*xA + rayVector[1]*yA + rayVector[2]*zA);
    double c = xA*xA + yA*yA + zA*zA - radius*radius;

    double disc  = b*b-4*a*c;
    if (disc < 0) // No imaginary answers
      return false;
    
    double sqrtDisc = sqrt(disc);
    double delta1   = (-b + sqrtDisc) / (2*a);
    double delta2   = (-b - sqrtDisc) / (2*a);

    if (delta1 >= 0) // First intersection valid
    {
      if ((delta2 < 0) || (delta1 < delta2)) // Second intersection invalid or farther away
        intersection = rayOrigin + delta1*rayVector;
      else // Second intersection valid and closer
        intersection = rayOrigin + delta2*rayVector;
      return true;
    }
    if (delta2 >= 0) // First intersection is invalid
      intersection = rayOrigin + delta2*rayVector;
    
    return false; // Both intersections are invalid
}

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

//----------------------------------------------------------------------------------------------------------------





namespace asp {
namespace isis {

class IsisInterfaceLineScanRot : public asp::isis::IsisInterfaceLineScan 
{
public:
  IsisInterfaceLineScanRot( std::string const& file ) : asp::isis::IsisInterfaceLineScan(file) {}
  virtual ~IsisInterfaceLineScanRot() {}

  /// Additional function to apply an in-camera rotation during this operation
  vw::Vector2
    point_to_pixel_rotated( vw::Vector3 const& point, vw::Vector3 const& rotAngles, int guessLine=-1) const;

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

  //std::cout << result << " ";
  return result;
}

Vector2
IsisInterfaceLineScanRot::point_to_pixel_rotated( vw::Vector3 const& point, vw::Vector3 const& rotAngles, int guessLine) const 
{

  // First seed LMA with an ephemeris time in the middle of the image
  double middle = lines() / 2;
  if (guessLine >= 0) // Use the input line as a guess if it was passed in
    middle = guessLine;
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
  if (status < 0)
  {
    //std::cout << "status = " << status << std::endl;
    //std::cout << "start_e   = " << start_e   << " solution = " << solution_e[0]        << std::endl;
    std::cout << "objective = " << objective << " result = "   << model(solution_e) << std::endl;
    std::cout << "Warning: Solver failed to find a solution!" << std::endl;
  }
  //VW_ASSERT( status >= 0,
  //           MathErr() << "Unable to project point into linescan camera , status = " << status );

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




#endif

