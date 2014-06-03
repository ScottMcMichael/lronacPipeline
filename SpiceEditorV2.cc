
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


/// \file SpiceEditor.cc
///

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp> 

#include <list>
#include <vector>
#include <string>
#include <iomanip>

//// CSpice include files
//#include "SpiceUsr.h"
//#include "SpiceZfc.h"

// ISIS files
#include <Cube.h>

#include <vw/Core/Exception.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Quaternion.h>

#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageView.h>

#include <stereo.h> // Using local version

#include <asp/IsisIO/IsisCameraModel.h>
#include <IsisInterfaceLineScanRot.h>


using namespace vw;
using std::endl;
using std::setprecision;
using std::setw;

static const int TRANSFORM_TYPE_GLOBAL   = 0; // Rotation only affects orientation
static const int TRANSFORM_TYPE_LOCAL    = 1; // Only apply a position offset in the spacecraft local frame
static const int TRANSFORM_TYPE_PC_ALIGN = 2; // Rotation affects position and orientation

struct Parameters : asp::BaseOptions 
{
  int          offsetCode;
  std::string  transformFile;
  int          transformType;
  
  bool debug;
  std::string sourceCubePath; // Used to obtain body position
};


bool handle_arguments(int argc, char* argv[], Parameters &opt) 
{ 
  std::string outputPrefix;
  po::options_description general_options("Options");
  general_options.add_options()
    ("debug",                  po::bool_switch(&opt.debug                 )->default_value(false),  "DEBUG mode")
    ("sourceCube",  po::value(&opt.sourceCubePath)->default_value(""), "Path to cube used to compute rotations")
    ("transformFile", po::value<std::string>(&opt.transformFile)->default_value(""), "Path to 3x4 matrix containing transform to apply (pc_align-style)")
    ("transformType", po::value<int>(&opt.transformType)->default_value(0),
          "Code to indicate how the transform is applied (0 = global, 1 = local(translation only), 2=pc_align");
    

  general_options.add( asp::BaseOptionsDescription(opt) );

  po::options_description positional("");
  positional.add_options();

  po::positional_options_description positional_desc;

  std::string usage("[options]\n");
  po::variables_map vm =
    asp::check_command_line( argc, argv, opt, general_options, general_options,
                             positional, positional_desc, usage );

  if (!vm.count("transformFile"))
    vw_throw( vw::ArgumentErr() << "Requires transform file in order to proceed.\n\n"
              << usage << general_options );
  
  return true;
}

/// Loads a transform file from pc_align into a rotation and translation matrix
bool loadThreeByFourTransform(const std::string &path, vw::Matrix3x3& R, vw::Vector3& T)
{
  // Open the file
  std::ifstream file(path.c_str());
  if (file.fail())
  {
    printf("Error reading transform file %s\n", path.c_str());
    return false;
  }

  // Read in one line at a time
  std::string line;
  for (int i=0; i<3; ++i)
  {
    std::getline(file, line);
    std::stringstream s(line); 
    s >> R[i][0] >> R[i][1] >> R[i][2] >> T[i]; // Fill in transform
  }

  // Clean up
  file.close();
  return true;
}


bool editSpiceFile(const Parameters &params)  
{
  //TODO: Figure out how to copy the cube to a new location!

  // First load the input transform data

  // Global transforms get stored here
  vw::Matrix3x3 correction_in_body_R;
  vw::Vector3   correction_in_body_T;
  // Local transforms get stored here
  vw::Matrix3x3 dummyRot;
  vw::Vector3 lronacOffset;

  if (params.transformType == TRANSFORM_TYPE_LOCAL)
  {
    if (!loadThreeByFourTransform(params.transformFile, dummyRot, lronacOffset))
      return false;
  }
  else // Global transform (from lronacAngleSolver or pc_align)
  {
    if (!loadThreeByFourTransform(params.transformFile, correction_in_body_R, correction_in_body_T))
      return false;
  }

  // Init an ISIS cube object
  //Isis::FileName cubePath(params.sourceCubePath.c_str());
  Isis::Cube cube;
  cube.open(params.sourceCubePath.c_str(), "rw");

  // Init an ASP interface to the cube also
  IsisInterfaceLineScanRot cubeInterface(params.sourceCubePath);

  // Make sure both information tables are there
  if (!cube.hasTable("InstrumentPosition")) {
    std::cout << "Error: Cube missing InstrumentPosition table!" << std::endl;
    return false;
  }
  if (!cube.hasTable("InstrumentPointing")) {
    std::cout << "Error: Cube missing InstrumentPointing table!" << std::endl;
    return false;
  }

  //----------------------------------------------------------------------------------------------
  //--- Rotation correction ---

  // If we are not applying a global transform there is nothing else to do here
  if (params.transformType != TRANSFORM_TYPE_LOCAL)
  {
    // Get the rotation table
    Isis::Table pointingTable("InstrumentPointing");
    cube.read(pointingTable);

    // Loop through all of the records
    int numRecords = pointingTable.Records();
    int numFields  = pointingTable.RecordFields(); // TODO: Verify the number of table fields
    for (int r=0; r<numRecords; ++r)
    {
      // Extract all desired parts of the record
      // - Fields: J2000Q0,J2000Q1,J2000Q2,J2000Q3,AV1,AV2,AV3,ET
      Isis::TableRecord& tableLine = pointingTable[r];
      double q0  = static_cast<double>(tableLine[0]); // Pointing quaternion in J2000 coordinates
      double q1  = static_cast<double>(tableLine[1]);
      double q2  = static_cast<double>(tableLine[2]);
      double q3  = static_cast<double>(tableLine[3]);
      double et  = static_cast<double>(tableLine[7]); // Ephemeris time
      
      //std::cout << " Old line values: ";
      //for (int c=0; c<numFields; ++c)
      //  std::cout << static_cast<double>(pointingTable[r][c]) << ", ";
      //std::cout << std::endl;
     
      // Convert input quaternion to rotation matrix
      // - Our correction is in body coordinates, not J2000 coordinates so we need to convert.
      vw::math::Quaternion<double> qIn(q0, q1, q2, q3);
      vw::Matrix3x3 instrument_from_J2000_R_in = qIn.rotation_matrix();

      // Get the body orientation
      vw::Matrix3x3 R_inst, body_from_J2000_R;
      cubeInterface.getMatricesAtTime(et, R_inst, body_from_J2000_R);

      // Go from J2000 coordinates to body coordinates
      vw::Matrix3x3 instrument_from_body_R = body_from_J2000_R * transpose(instrument_from_J2000_R_in);

      // Apply the input correction matrix
      vw::Matrix3x3 instrumentFixed_from_body_R = correction_in_body_R * instrument_from_body_R;

      // Convert back to J2000 frame
      vw::Matrix3x3 instrumentFixed_from_J2000 = transpose(transpose(body_from_J2000_R) * instrumentFixed_from_body_R);

      // Pack new rotation back into the cube via quaternion
      vw::math::Quaternion<double> qOut(instrumentFixed_from_J2000);
      tableLine[0] = qOut[0];
      tableLine[1] = qOut[1];
      tableLine[2] = qOut[2];
      tableLine[3] = qOut[3];
      pointingTable.Update(tableLine, r);
      
      //TODO: Update the table description

      //std::cout << " New line values: ";
      //for (int c=0; c<numFields; ++c)
      //{
      //  std::cout << static_cast<double>(pointingTable[r][c]) << ", ";
      //}
      //std::cout << std::endl;
    }
    printf("Writing rotation table!\n");
    cube.write(pointingTable);
  } // End of rotation correction in non-local transform case


  //----------------------------------------------------------------------------------------------
  //--- Position correction ---

  // Get the rotation table
  Isis::Table positionTable("InstrumentPosition");
  cube.read(positionTable);


  // Loop through all of the records
  int numRecords = positionTable.Records();
  int numFields  = positionTable.RecordFields(); // TODO: Verify the number of table fields
  for (int r=0; r<numRecords; ++r)
  {
    // Extract all desired parts of the record
    // - Fields: J2000X,J2000Y,J2000Z,J2000XV,J2000YV,J2000ZV,ET
    Isis::TableRecord& tableLine = positionTable[r];
    double x  = static_cast<double>(tableLine[0]); // Position in J2000 coordinates
    double y  = static_cast<double>(tableLine[1]);
    double z  = static_cast<double>(tableLine[2]);
    double et = static_cast<double>(tableLine[6]); // Ephemeris time
    
    //std::cout << " Old line values: ";
    //for (int c=0; c<numFields; ++c)
    //  std::cout << static_cast<double>(positionTable[r][c]) << ", ";
    //std::cout << std::endl;
   
    // The input coordinates
    vw::Vector3 instrument_from_J2000_T_in(x, y, z);

    // Get the instrument and body rotation at the same time
    vw::Matrix3x3 instrument_from_J2000_R, body_from_J2000_R;
    cubeInterface.getMatricesAtTime(et, instrument_from_J2000_R, body_from_J2000_R);

    vw::Vector3 newPosition;

    if (params.transformType == TRANSFORM_TYPE_LOCAL) // Apply existing lronac position offset
    {
      // Convert the LRONAC offset from meters to kilometers
      vw::Vector3 lronacOffsetKm = lronacOffset / 1000.0;
      
      // Now need to convert the LRONAC offset (spacecraft to frame) into J2000 coordinates
      // - Multiply by the inverted rotation matrix to go from spacecraft frame to J2000 frame
      vw::Vector3 instOffset_J2000Frame = transpose(instrument_from_J2000_R) * lronacOffsetKm;

      // Add the converted offset back in to the original coordinates (km to km here)
      newPosition = instrument_from_J2000_T_in + instOffset_J2000Frame;
    }
    else // Apply global transform from file
    {
      // Convert the body transform to the correct coordinate space
      vw::Matrix3x3 correction_in_J2000_R = correction_in_body_R * body_from_J2000_R;

      vw::Vector3   fixedBody_in_J2000_T;

      // Convert the state into units of meters
      vw::Vector3 stateMeters = instrument_from_J2000_T_in * 1000.0; // J2000 frame

      // Apply the corrective body-based transform to the vector currently represented in J2000 space
      vw::Vector3 stateBody = body_from_J2000_R * stateMeters; // Convert to body frame

      vw::Vector3 stateBodyFixed;
      if (params.transformType == TRANSFORM_TYPE_GLOBAL)
      {
        // In the transform that comes out of lronacAngleSolver (courtesy of the AdjustedCamera class),
        //  the solved for rotation does not affect the position at all so just add the translation.
        // -  AdjustedCameraModel applies transform about camera, not about the body!
        stateBodyFixed = stateBody + correction_in_body_T;
      }
      else // PC_ALIGN transform
      {
        // Apply the rotation to the position, then add in the translation
        stateBodyFixed = correction_in_body_R * stateBody;

        stateBodyFixed += correction_in_body_T;
      }

      // Now that the position is corrected in body frame, convert back to J2000 frame
      vw::Vector3 stateJ2000Fixed = transpose(body_from_J2000_R) * stateBodyFixed;
      newPosition = stateJ2000Fixed / 1000.0; // Convert from meters back to kilometers

    } // End of global transform case


    // Pack new position back into the cube
    tableLine[0] = newPosition[0];
    tableLine[1] = newPosition[1];
    tableLine[2] = newPosition[2];
    positionTable.Update(tableLine, r);
    
    //TODO: Update the table description

    //std::cout << " New line values: ";
    //for (int c=0; c<numFields; ++c)
    //{
    //  std::cout << static_cast<double>(positionTable[r][c]) << ", ";
    //}
    //std::cout << std::endl;
  }
  printf("Writing position table!\n");
  cube.write(positionTable);



  // Clean up
  cube.close();

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

    editSpiceFile(params);

  } ASP_STANDARD_CATCHES;

  printf("Exiting SpiceEditor program.\n");
  return 0;
}






