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


// ASP
#include <IsisInterfaceLineScanRot.h>
#include <iTime.h>
#include <vw/Math/EulerAngles.h>

using namespace vw;


// Construct
IsisInterfaceLineScanRot::IsisInterfaceLineScanRot(const std::string &filename ) : IsisInterface(filename), m_alphacube( *m_label ) {

  // Gutting Isis::Camera
  m_distortmap = m_camera->DistortionMap();
  m_focalmap   = m_camera->FocalPlaneMap();
  m_detectmap  = m_camera->DetectorMap();
}

// Custom Function to help avoid over invoking the deeply buried
// functions of Isis::Sensor
void IsisInterfaceLineScanRot::SetTime( Vector2 const& px, bool calc ) const {
  if ( px != m_c_location ) {
    m_c_location = px;
    m_detectmap->SetParent( m_alphacube.AlphaSample(px[0]),
                            m_alphacube.AlphaLine(px[1]) );

    if ( calc ) {
      // Calculating Spacecraft position and pose
      m_camera->instrumentPosition(&m_center[0]);
      m_center *= 1000;

      std::vector<double> rot_inst = m_camera->instrumentRotation()->Matrix();
      std::vector<double> rot_body = m_camera->bodyRotation()->Matrix();
      MatrixProxy<double,3,3> R_inst(&(rot_inst[0]));
      MatrixProxy<double,3,3> R_body(&(rot_body[0]));
      m_pose = Quat(R_body*transpose(R_inst));
    }
  }
}

/// Helper class for computing the projection line
class EphemerisLMA : public vw::math::LeastSquaresModelBase<EphemerisLMA> {
  vw::Vector3 m_point;
  Isis::Camera* m_camera;
  Isis::CameraDistortionMap *m_distortmap;
  Isis::CameraFocalPlaneMap *m_focalmap;
public:
  typedef vw::Vector<double> result_type; // Back project result
  typedef vw::Vector<double> domain_type; // Ephemeris time
  typedef vw::Matrix<double> jacobian_type;

  inline EphemerisLMA( vw::Vector3 const& point,
                       Isis::Camera* camera,
                       Isis::CameraDistortionMap* distortmap,
                       Isis::CameraFocalPlaneMap* focalmap ) 
                        : m_point(point), m_camera(camera), m_distortmap(distortmap), m_focalmap(focalmap) {}

  inline result_type operator()( domain_type const& x ) const;
};


// LMA for projecting point to linescan camera
EphemerisLMA::result_type
EphemerisLMA::operator()( EphemerisLMA::domain_type const& x ) const {

  // Setting Ephemeris Time
  m_camera->setTime(Isis::iTime(x[0]));

  // Calculating the look direction in camera frame
  Vector3 instru;
  m_camera->instrumentPosition(&instru[0]);
  instru *= 1000;  // Spice gives in km
  Vector3 lookB = normalize( m_point - instru );
  std::vector<double> lookB_copy(3);
  std::copy( lookB.begin(), lookB.end(), lookB_copy.begin() );
  std::vector<double> lookJ = m_camera->bodyRotation()->J2000Vector(lookB_copy);
  std::vector<double> lookC = m_camera->instrumentRotation()->ReferenceVector(lookJ);
  Vector3 look;
  std::copy( lookC.begin(), lookC.end(), look.begin() );

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
IsisInterfaceLineScanRot::point_to_pixel( Vector3 const& point ) const {

  // First seed LMA with an ephemeris time in the middle of the image
  double middle = lines() / 2;
  m_detectmap->SetParent( 1, m_alphacube.AlphaLine(middle) );
  double start_e = m_camera->time().Et();

  // Build LMA
  EphemerisLMA model( point, m_camera.get(), m_distortmap, m_focalmap );
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
  m_camera->instrumentPosition(&m_center[0]);
  m_center *= 1000;
  Vector3 look = normalize(point-m_center);

  // Calculating Rotation to camera frame
  std::vector<double> rot_inst = m_camera->instrumentRotation()->Matrix();
  std::vector<double> rot_body = m_camera->bodyRotation()->Matrix();
  MatrixProxy<double,3,3> R_inst(&(rot_inst[0]));
  MatrixProxy<double,3,3> R_body(&(rot_body[0]));
  m_pose = Quat(R_body*transpose(R_inst));

  look = inverse(m_pose).rotate( look );
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
                          //m_rotation(rotation[0], rotation[1], 0) // Only use certain rotations
                          {
                          }

  inline result_type operator()( domain_type const& x ) const;
};


// LMA for projecting point to linescan camera
EphemerisLMA_rot::result_type
EphemerisLMA_rot::operator()( EphemerisLMA::domain_type const& x ) const {

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



/// Hack function to insert an additional camera rotation into this function 
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
  Vector3 offsetAngles =  rotAngles;
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

bool 
IsisInterfaceLineScanRot::getMatricesAtTime(const double et, Matrix3x3 &R_inst, Matrix3x3 &R_body)
{
  // Converting now to pixel
  m_camera->setTime(Isis::iTime(et));

  // Calculating Rotation to camera frame
  std::vector<double> rot_inst = m_camera->instrumentRotation()->Matrix();
  std::vector<double> rot_body = m_camera->bodyRotation()->Matrix();
  R_inst = Matrix3x3(&(rot_inst[0])); // Rotation of spacecraft and instrument relative to J2000
  R_body = Matrix3x3(&(rot_body[0])); // Rotation of planet relative to J2000 (earth-centered) frame
  return true;
} 
  


Vector3
IsisInterfaceLineScanRot::pixel_to_vector( Vector2 const& pix ) const {
  Vector2 px = pix + Vector2(1,1);
  SetTime( px, true );

  // Projecting to get look direction
  Vector3 result;
  m_focalmap->SetDetector( m_detectmap->DetectorSample(),
                           m_detectmap->DetectorLine() );
  m_distortmap->SetFocalPlane( m_focalmap->FocalPlaneX(),
                               m_focalmap->FocalPlaneY() );
  result[0] = m_distortmap->UndistortedFocalPlaneX();
  result[1] = m_distortmap->UndistortedFocalPlaneY();
  result[2] = m_distortmap->UndistortedFocalPlaneZ();
  result = normalize( result );
  result = m_pose.rotate(result);
  return result;
}

Vector3
IsisInterfaceLineScanRot::camera_center( Vector2 const& pix ) const {
  Vector2 px = pix + Vector2(1,1);
  SetTime( px, true );
  return m_center;
}

Quat
IsisInterfaceLineScanRot::camera_pose( Vector2 const& pix ) const {
  Vector2 px = pix + Vector2(1,1);
  SetTime( px, true );
  return m_pose;
}
