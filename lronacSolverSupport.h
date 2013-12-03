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

//###############################################################################################
//###############################################################################################
//###############################################################################################
// Ray/Sphere intersection code

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

//###############################################################################################
//###############################################################################################
//###############################################################################################
// IsisInterfaceLineScanRot - Class replacment to insert local rotation

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




//###############################################################################################
//###############################################################################################
//###############################################################################################
// IsisInterfaceLineScanRot - Class replacment to insert local rotation

namespace vw {
namespace camera {

/// This class is useful if you have an existing camera model, and
/// you want to systematically "tweak" its extrinsic parameters
/// (position and pose).  This is particularly useful in Bundle
/// Adjustment.
class AdjustedCameraModelRot : public CameraModel 
{
private:
  
  boost::shared_ptr<IsisInterfaceLineScanRot> m_camera;
  Vector3 m_translation;
  Quat m_rotation;
  Quat m_rotation_inverse;

public:
  AdjustedCameraModelRot(){}
  
  AdjustedCameraModelRot(boost::shared_ptr<IsisInterfaceLineScanRot> camera_model) : m_camera(camera_model) 
  {
    m_rotation = Quat(math::identity_matrix<3>());
    m_rotation_inverse = Quat(math::identity_matrix<3>());
  }

  AdjustedCameraModelRot(boost::shared_ptr<IsisInterfaceLineScanRot> camera_model,
                      Vector3 const& translation, Quat const& rotation) :
    m_camera(camera_model), m_translation(translation), m_rotation(rotation), m_rotation_inverse(inverse(rotation)) {}

  virtual ~AdjustedCameraModelRot() {}
  virtual std::string type() const { return "Adjusted"; }

  Vector3 translation() const { return m_translation; }
  Quat rotation() const { return m_rotation; }
  Matrix<double,3,3> rotation_matrix() const { return m_rotation.rotation_matrix(); }
  //Vector3 axis_angle_rotation() const;
  //void set_rotation(Quat const&);

  template <class MatrixT>
  void set_rotation(MatrixBase<MatrixT> const& m) 
  {
    m_rotation = Quat(m.impl());
    m_rotation_inverse = inverse(m_rotation);
  }
  template <class VectorT>
  void set_translation(VectorBase<VectorT> const& v) 
  {
    m_translation = v.impl();
  }
  template <class VectorT>
  void set_axis_angle_rotation(VectorBase<VectorT> const& v) 
  {
    this->set_rotation( axis_angle_to_quaternion(v.impl()) );
  }

  //virtual Vector2 point_to_pixel (Vector3 const&) const;
  //virtual Vector3 pixel_to_vector (Vector2 const&) const;
  //virtual Vector3 camera_center (Vector2 const&) const;
  //virtual Quat camera_pose(Vector2 const&) const;

  //void write(std::string const&);
  //void read(std::string const&);

  //friend std::ostream& operator<<(std::ostream&, AdjustedCameraModel const&);
  
Vector3 axis_angle_rotation() const {
  Quat quat = this->rotation();
  return quat.axis_angle();
}

void set_rotation(Quat const& rotation) {
  m_rotation = rotation;
  m_rotation_inverse = inverse(m_rotation);
}

Vector2 point_to_pixel (Vector3 const& point) const {
  Vector3 offset_pt = point-m_camera->camera_center(Vector2(0,0))-m_translation;
  Vector3 new_pt = m_rotation_inverse.rotate(offset_pt) + m_camera->camera_center(Vector2(0,0));
  return m_camera->point_to_pixel(new_pt);
}

/// New code to support passing local rotation angles into this function
Vector2 point_to_pixel_rotated( vw::Vector3 const& point, vw::Vector3 const& rotAngles, int guessLine) const
{
  Vector3 offset_pt = point-m_camera->camera_center(Vector2(0,0))-m_translation;
  Vector3 new_pt = m_rotation_inverse.rotate(offset_pt) + m_camera->camera_center(Vector2(0,0));
  return m_camera->point_to_pixel_rotated(new_pt, rotAngles, guessLine);
}

Vector3 pixel_to_vector (Vector2 const& pix) const {
  return m_rotation.rotate(m_camera->pixel_to_vector(pix));
}

Vector3 camera_center(Vector2 const& pix) const {
  return m_camera->camera_center(pix) + m_translation;
}

Quat camera_pose(Vector2 const& pix) const {
  return m_rotation*m_camera->camera_pose(pix);
}
/*
void write(std::string const& filename) {
  std::ofstream ostr(filename.c_str());
  ostr << m_translation[0] << " " << m_translation[1]
       << " " << m_translation[2] << "\n";
  ostr << m_rotation.w() << " " << m_rotation.x() << " "
       << m_rotation.y() << " " << m_rotation.z() << "\n";
}

void read(std::string const& filename) {
  Vector4 c;
  Vector3 pos;
  std::ifstream istr(filename.c_str());
  istr >> pos[0] >> pos[1] >> pos[2];
  istr >> c[0] >> c[1] >> c[2] >> c[3];
  this->set_translation(pos);
  this->set_rotation(Quat(c));
}
*/
/*
std::ostream& camera::operator<<(std::ostream& ostr,
           AdjustedCameraModel const& cam ) {
  ostr << "AdjustedCameraModel(Trans: " << cam.m_translation << " Rot: "
       << cam.m_rotation << " Cam: " << cam.m_camera->type() << ")\n";
  return ostr;
}
*/
  
};
  
} // end namespace camera
} // end namespace vw

#endif

