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
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Stereo/PreFilter.h>
#include <vw/Stereo/Correlate.h>
#include <vw/Stereo/StereoModel.h>
#include <asp/IsisIO/IsisCameraModel.h> //::point_to_pixel>
#include <vw/Math/LevenbergMarquardt.h>
#include <vw/Math/EulerAngles.h>

#include <asp/Core/IntegralAutoGainDetector.h>
//#include <asp/IsisIO/IsisInterfaceLineScan.h>
#include <IsisInterfaceLineScanRot.h>

#include <stereo.h>

//###############################################################################################
//###############################################################################################
//###############################################################################################
// Ray/Sphere intersection code

//TODO: Replace this code!!!!


bool raySphereIntersect(const vw::Vector3 &rayOrigin, const vw::Vector3 &rayVector, const double radius, vw::Vector3 &intersection)
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


/// Expands the original AdjustedCameraModel class to allow global and local rotations
class AdjustedCameraModelRot : public LocalRotCameraModel
{
private:
  
  boost::shared_ptr<LocalRotCameraModel> m_camera;
  vw::Vector3 m_translation;
  vw::Quat    m_rotation;
  vw::Quat    m_rotation_inverse;

public:
  AdjustedCameraModelRot(){}
  
  AdjustedCameraModelRot(boost::shared_ptr<LocalRotCameraModel> camera_model) : m_camera(camera_model)
  {
    m_rotation         = vw::Quat(vw::math::identity_matrix<3>());
    m_rotation_inverse = vw::Quat(vw::math::identity_matrix<3>());
  }

  AdjustedCameraModelRot(boost::shared_ptr<LocalRotCameraModel> camera_model,
                      vw::Vector3 const& translation, vw::Quat const& rotation) :
    m_camera(camera_model), m_translation(translation), m_rotation(rotation), m_rotation_inverse(inverse(rotation)) {}

  virtual ~AdjustedCameraModelRot() {}
  virtual std::string type() const { return "Adjusted"; }

  vw::Vector3            translation    () const { return m_translation; }
  vw::Quat               rotation       () const { return m_rotation; }
  vw::Matrix<double,3,3> rotation_matrix() const { return m_rotation.rotation_matrix(); }

  template <class MatrixT>
  void set_rotation(vw::MatrixBase<MatrixT> const& m)
  {
    m_rotation = Quat(m.impl());
    m_rotation_inverse = inverse(m_rotation);
  }
  template <class VectorT>
  void set_translation(vw::VectorBase<VectorT> const& v)
  {
    m_translation = v.impl();
  }
  template <class VectorT>
  void set_axis_angle_rotation(vw::VectorBase<VectorT> const& v)
  {
    this->set_rotation( axis_angle_to_quaternion(v.impl()) );
  }
  
  vw::Vector3 axis_angle_rotation() const {
  vw::Quat quat = this->rotation();
  return quat.axis_angle();
}

void set_rotation(vw::Quat const& rotation) {
  m_rotation = rotation;
  m_rotation_inverse = inverse(m_rotation);
}

virtual vw::Vector2 point_to_pixel (vw::Vector3 const& point) const {
  vw::Vector3 offset_pt = point-m_camera->camera_center(vw::Vector2(0,0))-m_translation;
  vw::Vector3 new_pt    = m_rotation_inverse.rotate(offset_pt) + m_camera->camera_center(vw::Vector2(0,0));
  return m_camera->point_to_pixel(new_pt);
}

/// New code to support passing local rotation angles into this function
virtual vw::Vector2 point_to_pixel_rotated( vw::Vector3 const& point, vw::Vector3 const& rotAngles, int guessLine) const
{
  vw::Vector3 offset_pt = point-m_camera->camera_center(vw::Vector2(0,0))-m_translation;
  vw::Vector3 new_pt    = m_rotation_inverse.rotate(offset_pt) + m_camera->camera_center(vw::Vector2(0,0));
  return m_camera->point_to_pixel_rotated(new_pt, rotAngles, guessLine);
}

virtual vw::Vector3 pixel_to_vector (vw::Vector2 const& pix) const {
  return m_rotation.rotate(m_camera->pixel_to_vector(pix));
}

/// New code to support passing local rotation angles into this function
virtual vw::Vector3 pixel_to_vector_rotated (vw::Vector2 const& pix, vw::Vector3 const& rotAngles) const {
  return m_rotation.rotate(m_camera->pixel_to_vector_rotated(pix, rotAngles));
}

virtual vw::Vector3 camera_center(vw::Vector2 const& pix) const {
  return m_camera->camera_center(pix) + m_translation;
}

virtual vw::Quat camera_pose(vw::Vector2 const& pix) const {
  return m_rotation*m_camera->camera_pose(pix);
}
  
};
  
#endif

