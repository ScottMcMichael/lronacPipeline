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


/// \file IsisInterfaceLineScane.h
///
/// Line Scan Camera Interface with ISIS
///
#ifndef __LRONACPIPELINE_INTERFACE_LINESCAN_H__
#define __LRONACPIPELINE_INTERFACE_LINESCAN_H__

// VW & ASP
#include <vw/Math/LevenbergMarquardt.h>
#include <vw/Camera/CameraModel.h>
#include <asp/IsisIO/IsisInterface.h>

// Isis
#include <CameraDetectorMap.h>
#include <CameraDistortionMap.h>
#include <CameraFocalPlaneMap.h>
#include <AlphaCube.h>


  class IsisInterfaceLineScanRot : public asp::isis::IsisInterface, public vw::camera::CameraModel {

  public:
    IsisInterfaceLineScanRot(const std::string &file );

    virtual ~IsisInterfaceLineScanRot() {}

    virtual std::string type()  { return "LineScan"; }
    virtual std::string type() const { return "LineScan"; }

    // Standard Methods
    //-------------------------------------------------

    virtual vw::Vector2
      point_to_pixel( vw::Vector3 const& point ) const;
    virtual vw::Vector3
      pixel_to_vector( vw::Vector2 const& pix ) const;
    virtual vw::Vector3
      camera_center( vw::Vector2 const& pix = vw::Vector2(1,1) ) const;
    virtual vw::Quat
      camera_pose( vw::Vector2 const& pix = vw::Vector2(1,1) ) const;

    /// Additional function to apply an in-camera rotation during this operation
    vw::Vector2
      point_to_pixel_rotated( vw::Vector3 const& point, vw::Vector3 const& rotAngles, int guessLine=-1) const;

    /// Returns the instrument and body matrices at the given time
    bool getMatricesAtTime(const double et, vw::Matrix3x3 &R_inst, vw::Matrix3x3 &R_body);

  protected:

    // Custom Variables
    Isis::CameraDistortionMap *m_distortmap;
    Isis::CameraFocalPlaneMap *m_focalmap;
    Isis::CameraDetectorMap   *m_detectmap;
    mutable Isis::AlphaCube   m_alphacube; // Doesn't use const

  private:

    // Custom Fuctions
    mutable vw::Vector2 m_c_location;
    mutable vw::Vector3 m_center;
    mutable vw::Quat m_pose;
    void SetTime( vw::Vector2 const& px,
                  bool calc=false ) const;
  };



#endif//__LRONACPIPELINE_INTERFACE_LINESCAN_H__
