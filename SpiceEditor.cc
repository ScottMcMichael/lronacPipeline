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

// CSpice include files
#include "SpiceUsr.h"
#include "SpiceZfc.h"

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


struct Parameters : asp::BaseOptions 
{
  int          offsetCode;
  std::string  spkDataOutputPath;
  std::string  ckDataOutputPath;
  std::vector<std::string>  kernelPaths;
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
    ("outputPrefix",  po::value(&outputPrefix)->default_value(""), "Output prefix")
    ("kernels", po::value<std::vector<std::string> >(&opt.kernelPaths)->multitoken(), "Paths to all required kernel files")
    ("transformFile", po::value<std::string>(&opt.transformFile)->default_value(""), "Path to 3x4 matrix containing transform to apply (pc_align-style)")
    ("transformType", po::value<int>(&opt.transformType)->default_value(0), "Code to indicate how the transform is applied (0 = global, 1 = local = translation only!)");
    

  general_options.add( asp::BaseOptionsDescription(opt) );

  po::options_description positional("");
  positional.add_options();

  po::positional_options_description positional_desc;

  std::string usage("[options] ");
  po::variables_map vm =
    asp::check_command_line( argc, argv, opt, general_options, general_options,
                             positional, positional_desc, usage );

  if ( !vm.count("kernels") || !vm.count("outputPrefix") || !vm.count("transformFile"))
    vw_throw( vw::ArgumentErr() << "Requires kernels, transform file, and output path in order to proceed.\n\n"
              << usage << general_options );


  // Finish paths to CK and SPK output data files
  opt.spkDataOutputPath = outputPrefix + "-spkData.txt";
  opt.ckDataOutputPath  = outputPrefix + "-ckData.txt";
  
  return true;
}

/// Loads a transform file from pc_align into a rotation and translation matrix
bool loadThreeByFourTransform(const std::string &path, SpiceDouble R[][3], SpiceDouble *T)
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
  const SpiceInt    MOON_CODE           = 301;
  const SpiceInt    LRO_CLOCK_ID        = -85;
  const SpiceInt    LRO_SPACECRAFT_CODE = -85000; // Instrument ID code
  const std::string J2000_FRAME_STRING  = "J2000";
  const std::string MOON_FRAME_STRING   = "IAU_MOON";
  
  const std::string MOON_STRING         = "moon";
  const std::string EARTH_STRING        = "earth";

  // Output kernel files have this many entries
  // - Could definately do something better than hard-coding this number
  const int NUM_STEPS = 5000;
  
  //TODO: Read these from file
  // Target and observer are common for all operations
  std::string target   = MOON_STRING;
  std::string observer = EARTH_STRING;
  std::string absCorr  = "NONE";

  // Find the ck and spk files TODO: Handle multiple!
  std::string spkFile = "";
  std::string ckFile  = "";
  for (size_t i=0; i<params.kernelPaths.size(); ++i)
  {
    // Check correct extension and not in target folder
    if ( (params.kernelPaths[i].find("/tspk/") == std::string::npos) &&
	 (params.kernelPaths[i].find(".bsp") != std::string::npos)    )
    {
      spkFile = params.kernelPaths[i];
      printf("Found spk file %s\n", spkFile.c_str());
    }
    // Check correct extension
    if ( /*(params.kernelPaths[i].find("/ck/") != std::string::npos) &&*/
	 (params.kernelPaths[i].find(".bc") != std::string::npos)    )
    {
      ckFile = params.kernelPaths[i];
      printf("Found ck file %s\n", ckFile.c_str());
    }  
  }
  if (spkFile.empty())
  {
    printf("Failed to find spk file!\n");
    return false;
  }  
  if (ckFile.empty())
  {
    printf("Failed to find ck file!\n");
    return false;
  }
  
  // Load all of the kernel files
  printf("Loading all source kernels (this can take a while)\n");
  for (size_t i=0; i<params.kernelPaths.size(); ++i)
    furnsh_c(params.kernelPaths[i].c_str());

  
  // Try loading the source cube
  double minSourceCubeEt=0, maxSourceCubeEt=-1;
  boost::shared_ptr<IsisInterfaceLineScanRot> sourceCubePtr;
  if (!params.sourceCubePath.empty())
  {
    printf("Loading source cube file %s\n", params.sourceCubePath.c_str());
    sourceCubePtr.reset(new IsisInterfaceLineScanRot(params.sourceCubePath));
    int    numLines   = sourceCubePtr->lines();
    int    numSamples = sourceCubePtr->samples();
    minSourceCubeEt = sourceCubePtr->ephemeris_time(vw::Vector2(numSamples/2, 0));
    maxSourceCubeEt = sourceCubePtr->ephemeris_time(vw::Vector2(numSamples/2, numLines-1));
  }

  // Load the transform to apply  
  const bool localTransform = (params.transformType > 0);

  // Global transforms get stored here
  SpiceDouble planetFixed_from_planet_R[3][3];
  SpiceDouble planetFixed_from_planet_T[3];
  if (!localTransform)
  {
    if (!loadThreeByFourTransform(params.transformFile, planetFixed_from_planet_R, planetFixed_from_planet_T))
      return false;
     
  }
  
  // Local transforms get stored here
  SpiceDouble dummyRot[3][3];
  SpiceDouble lronacOffset[3];
  if (localTransform)
  {
    if (!loadThreeByFourTransform(params.transformFile, dummyRot, lronacOffset))
      return false;
  }


  if (params.debug) // Print out a test case then quit
  {
    SpiceDouble et = 342594244.951690;
    printf("et = %lf\n", et);

    SpiceDouble  tol = 0; // Make sure we can do this 

    // Convert ephemeris time to spacecraft clock time
    SpiceDouble sclkdp;
    sce2c_c(LRO_CLOCK_ID, et, &sclkdp);

    // Try to get the spacecraft orientation at that time (spacecraft_from_J2000)
    SpiceDouble  spacecraft_from_J2000_R[3][3];
    SpiceDouble  clkout;
    SpiceBoolean orientationFound;
    ckgp_c(LRO_SPACECRAFT_CODE, sclkdp, tol, J2000_FRAME_STRING.c_str(), spacecraft_from_J2000_R, &clkout, &orientationFound);
    if (!orientationFound)
    {
      printf("SpiceEditor Error: Failed to obtain spacecraft orientation data!!!!!!!!!!!!!!!!\n");
      return false;
    }

    // Try to get the planet orientation at that time (planet_from_J2000)
    // - [this matrix] * [j2000 vector] = [moon vector]
    SpiceDouble  planet_from_J2000_R[3][3];

    if ((et >= minSourceCubeEt) && (et <= maxSourceCubeEt)) // If this time falls within the cube time
    {
      // Get the planet orientation from the source cube
      vw::Matrix3x3 R_inst, R_body;
      sourceCubePtr->getMatricesAtTime(et, R_inst, R_body);
      for (int r=0; r<3; ++r)
      {
        for (int c=0; c<3; ++c)
        {
          planet_from_J2000_R[r][c] = R_body[r][c];
        }
      }
    }
    else // Get the planet orientation from NAIF calls
    {
      pxform_c(J2000_FRAME_STRING.c_str(), MOON_FRAME_STRING.c_str(), et, planet_from_J2000_R);
//    tipbod_c(J2000_FRAME_STRING.c_str(), MOON_CODE, et, planet_from_J2000_R); // Seems equivalent
    }


//    SpiceDouble  J2000_from_planet_R[3][3];
//    pxform_c(MOON_FRAME_STRING.c_str(), J2000_FRAME_STRING.c_str(), et, J2000_from_planet_R);

    // Convert the planet transform to the correct coordinate space
    SpiceDouble planetFixed_from_J2000_R[3][3];

    // Have fixedPlanet_from_planet (or reverse), sc and planet from J2000

//        // Convert correcting rotation and translation into J2000 frame from planet frame
//        mxm_c(planetFixed_from_planet_R, planet_from_J2000_R, planetFixed_from_J2000_R);

    // Apply the transform
    SpiceDouble spacecraftFixed_from_J2000_R[3][3];
    SpiceDouble spacecraft_from_Planet_R[3][3];
    SpiceDouble spacecraftFixed_from_Planet_R[3][3];

    // Is this actually reversed?
    mxmt_c(planet_from_J2000_R,      spacecraft_from_J2000_R,      spacecraft_from_Planet_R);       // Convert to moon frame
    //mxmt_c(spacecraft_from_J2000_R,      planet_from_J2000_R,      spacecraft_from_Planet_R);       // Convert to moon frame

/*
    SpiceDouble spacecraft_from_Planet_R_2[3][3];
    pxform_c(MOON_FRAME_STRING.c_str(), "LRO_SC_BUS", et, spacecraft_from_Planet_R_2);
    printf("\nConverted to planet 2!!! (SC in planet)\n");
    for (int q=0; q<3; ++q)
    {
      for (int s=0; s<3; ++s)
        std::cout << " " << spacecraft_from_Planet_R_2[q][s]; // Write out the new rotation
      std::cout << std::endl;
    }
*/
    
    mxm_c(planetFixed_from_planet_R,     spacecraft_from_Planet_R, spacecraftFixed_from_Planet_R);  // Apply correction
    
    SpiceDouble temp_R[3][3];
    mtxm_c(planet_from_J2000_R,      spacecraftFixed_from_Planet_R,      temp_R);
    xpose_c(temp_R, spacecraftFixed_from_J2000_R);
    //mxm_c(J2000_from_planet_R,      spacecraftFixed_from_Planet_R, spacecraftFixed_from_J2000_R); // Return to J2000 frame
    //mxm_c(spacecraftFixed_from_Planet_R, planet_from_J2000_R,      spacecraftFixed_from_J2000_R); // Return to J2000 frame

    printf("\nOld rotation matrix (SC in J2000)\n");
    for (int q=0; q<3; ++q)
    {
      for (int s=0; s<3; ++s)
        std::cout << " " << spacecraft_from_J2000_R[q][s]; // Write out the new rotation
      std::cout << std::endl;
    }

    printf("\nConversion rotation matrix (planet in J2000)\n");
    for (int q=0; q<3; ++q)
    {
      for (int s=0; s<3; ++s)
        std::cout << " " << planet_from_J2000_R[q][s]; // Write out the new rotation
      std::cout << std::endl;
    }

    printf("\nConverted to planet (SC in planet)\n");
    for (int q=0; q<3; ++q)
    {
      for (int s=0; s<3; ++s)
        std::cout << " " << spacecraft_from_Planet_R[q][s]; // Write out the new rotation
      std::cout << std::endl;
    }

    printf("\nFixed in planet (new SC in planet)\n");
    for (int q=0; q<3; ++q)
    {
      for (int s=0; s<3; ++s)
        std::cout << " " << spacecraftFixed_from_Planet_R[q][s]; // Write out the new rotation
      std::cout << std::endl;
    }


    printf("\nModifying rotation matrix\n");
    for (int q=0; q<3; ++q)
    {
      for (int s=0; s<3; ++s)
        std::cout << " " << planetFixed_from_planet_R[q][s]; // Write out the new rotation
      std::cout << std::endl;
    }     

    printf("\nNew rotation matrix (new SC in J2000)\n");
    for (int q=0; q<3; ++q)
    {
      for (int s=0; s<3; ++s)
        std::cout << " " << spacecraftFixed_from_J2000_R[q][s]; // Write out the new rotation
      std::cout << std::endl;
    }

    return true;
  }

  SpiceChar timeString[51];
  SPICEDOUBLE_CELL(cover, 2000);
  SpiceDouble b, e;  
  
  // Check out the coverage of the ck (orientation) file
  //printf("CK file coverage:\n");
  SPICEINT_CELL(ckIds, 1000);
  ckobj_c (ckFile.c_str(), &ckIds); // Get list of CK id's in the file
  
  //TODO: Only operate on the CK ID of interest!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  for (int i=0; i<card_c(&ckIds); ++i) // For each ID in the file
  {
    SpiceInt bodyId = SPICE_CELL_ELEM_I(&ckIds, i);
    
    // Get the coverage window
    SpiceBoolean  needav = 0;
    scard_c(0, &cover);
    ckcov_c(ckFile.c_str(), bodyId, needav, "SEGMENT", 0.0, "TDB", &cover);
       
   
    // Get the number of intervals in the coverage window.
    SpiceInt numIntervals = wncard_c(&cover);    

    //printf ( "\nCoverage for CK object %ld\n", (long int)bodyId );
    
    // Convert the coverage interval start and stop times to TDB calendar strings.
    for (int j=0; j<numIntervals; j++)
    {
      // Get the endpoints of the jth interval.
      wnfetd_c(&cover, j, &b, &e);

      // Convert the endpoints to TDB calendar format time strings and display them.
      timout_c(b, "YYYY MON DD HR:MN:SC.### (TDB) ::TDB", 51, timeString);

      //printf("\nInterval:  %ld\nStart:     %s\n", j, timeString);

      //timout_c(e, "YYYY MON DD HR:MN:SC.### (TDB) ::TDB", 51, timeString);
      //printf ( "Stop:      %s\n", timeString );


      //printf("b = %lf, e = %lf\n", b, e); 

      // If we are not applying a global transform there is nothing else to do here
      if (localTransform)
        continue;
        
      // Correct the rotation at even intervals
      SpiceDouble startEt = b+2; // Hack to avoid different start times
      SpiceDouble stopEt  = e-2;

      if (maxSourceCubeEt > minSourceCubeEt) // Limit processing to cube duration
      {
        startEt = minSourceCubeEt - 2;
        stopEt  = maxSourceCubeEt + 2;
      }
      SpiceDouble stepSize = (stopEt - startEt) / NUM_STEPS;
    
      std::ofstream outputFile;
      printf("Writing file %s\n", params.ckDataOutputPath.c_str());
      outputFile.open(params.ckDataOutputPath.c_str());

      // For the number of specified steps
      SpiceDouble lastEt;
      for (int i=0; i<NUM_STEPS; ++i)
      {

        SpiceDouble state[6];
        SpiceDouble lightTime;
        SpiceDouble et = startEt + stepSize*(SpiceDouble)i;
        lastEt = et;

        SpiceDouble  tol = 0; // Make sure we can do this 

        // Convert ephemeris time to spacecraft clock time
        SpiceDouble sclkdp;
        sce2c_c(LRO_CLOCK_ID, et, &sclkdp);

        // Try to get the spacecraft orientation at that time (spacecraft_from_J2000)
        SpiceDouble  spacecraft_from_J2000_R[3][3];
        SpiceDouble  clkout;
        SpiceBoolean orientationFound;
        ckgp_c(LRO_SPACECRAFT_CODE, sclkdp, tol, J2000_FRAME_STRING.c_str(), spacecraft_from_J2000_R, &clkout, &orientationFound);
        if (!orientationFound)
        {
          printf("SpiceEditor Error: Failed to obtain spacecraft orientation data at et %lf!!!!!!!!!!!!!!!!\n", et);
          return false;
        }

        // Try to get the planet orientation at that time (planet_from_J2000)
        // - [this matrix] * [j2000 vector] = [moon vector]
        SpiceDouble  planet_from_J2000_R[3][3];
        pxform_c(J2000_FRAME_STRING.c_str(), MOON_FRAME_STRING.c_str(), et, planet_from_J2000_R);

        SpiceDouble  J2000_from_planet_R[3][3];
        if ((et >= minSourceCubeEt) && (et <= maxSourceCubeEt)) // If this time falls within the cube time
        {
          // Get the planet orientation from the source cube
          // - For some reason this produces different results than going straight to the NAIF functions
          vw::Matrix3x3 R_inst, R_body;
          sourceCubePtr->getMatricesAtTime(et, R_inst, R_body);
          for (int r=0; r<3; ++r)
          {
            for (int c=0; c<3; ++c)
            {
              planet_from_J2000_R[r][c] = R_body[r][c];
            }
          }
        }
        else // Get the planet orientation from NAIF calls
        {
          pxform_c(J2000_FRAME_STRING.c_str(), MOON_FRAME_STRING.c_str(), et, planet_from_J2000_R);
        }


        // Convert the planet transform to the correct coordinate space
        SpiceDouble planetFixed_from_J2000_R[3][3];

        // Have fixedPlanet_from_planet (or reverse), sc and planet from J2000

//        // Convert correcting rotation and translation into J2000 frame from planet frame
//        mxm_c(planetFixed_from_planet_R, planet_from_J2000_R, planetFixed_from_J2000_R);

        // Apply the transform
        SpiceDouble spacecraftFixed_from_J2000_R [3][3];
        SpiceDouble spacecraft_from_Planet_R     [3][3];
        SpiceDouble spacecraftFixed_from_Planet_R[3][3];
        
        
        mxmt_c(planet_from_J2000_R,      spacecraft_from_J2000_R,      spacecraft_from_Planet_R);       // Convert to moon frame
        
        mxm_c(planetFixed_from_planet_R,  spacecraft_from_Planet_R,     spacecraftFixed_from_Planet_R);  // Apply correction
        
        // Convert back to J2000 frame
        SpiceDouble temp_R[3][3];
        mtxm_c(planet_from_J2000_R,      spacecraftFixed_from_Planet_R,  temp_R);
        xpose_c(temp_R, spacecraftFixed_from_J2000_R);

        

        // Dump ET, new rotation... line to a text file
        outputFile.precision(16);
        outputFile << et;
        for (int q=0; q<3; ++q)
          for (int s=0; s<3; ++s)
            outputFile << " " << spacecraftFixed_from_J2000_R[q][s]; // Write out the new rotation
        outputFile << std::endl;
      } // End of steps loop

      outputFile.close(); // Close read data file

    } // End of loop through intervals

  } // End of loop through bodies

  //printf("SPK file coverage:\n");

  // Get list of bodies covered by the kernel
  SPICEINT_CELL(ids, 1000);
  spkobj_c(spkFile.c_str(), &ids);


  //TODO: Just loop through the times we pull from the tabledump

  // For each body
  for (int i=0; i<card_c(&ids); ++i)
  {
    SpiceInt bodyId = SPICE_CELL_ELEM_I(&ids, i);

    // Get the coverage window
    scard_c(0, &cover);
    spkcov_c(spkFile.c_str(), bodyId, &cover);


    // Get the number of intervals in the coverage window.
    SpiceInt numIntervals = wncard_c(&cover);    

    //printf("\nCoverage for SPK object %ld\n", (long int)bodyId);

    // Convert the coverage interval start and stop times to TDB calendar strings.
    for (int j=0; j<numIntervals; j++)
    {
      // Get the endpoints of the jth interval.
      wnfetd_c(&cover, j, &b, &e);

      // Convert the endpoints to TDB calendar format time strings and display them.
      timout_c(b, "YYYY MON DD HR:MN:SC.### (TDB) ::TDB", 51, timeString);

      //printf("\nInterval:  %ld\nStart:     %s\n", j, timeString);

      timout_c(e, "YYYY MON DD HR:MN:SC.### (TDB) ::TDB", 51, timeString);
      //printf("Stop:      %s\n", timeString);

      //printf("b = %lf, e = %lf\n", b, e); 
      // Modify the position at even intervals
      SpiceDouble startEt = b+2; // Hack to avoid different start times
      SpiceDouble stopEt  = e-2;

      if (maxSourceCubeEt > minSourceCubeEt) // Limit processing to cube duration
      {
        startEt = minSourceCubeEt - 2;
        stopEt  = maxSourceCubeEt + 2;
      }
      SpiceDouble stepSize = (stopEt - startEt) / NUM_STEPS;

      std::ofstream outputFile;
      printf("Writing file %s\n", params.spkDataOutputPath.c_str());
      outputFile.open(params.spkDataOutputPath.c_str());

      // For the number of specified steps
      SpiceDouble lastEt;
      for (int i=0; i<NUM_STEPS; ++i)
      {

        SpiceDouble state[6];
        SpiceDouble lightTime;
        SpiceDouble et = startEt + stepSize*(SpiceDouble)i;
        lastEt = et;

        // Retrieve the position of the spacecraft at this time
        spkez_c(bodyId,  et,  J2000_FRAME_STRING.c_str(), absCorr.c_str(), MOON_CODE, state, &lightTime); // LRO relative to Moon in J2000 frame, units are kilometers and km/sec

        SpiceDouble tol = 0; // Make sure we can do this 

        // Convert ephemeris time to spacecraft clock time
        SpiceDouble sclkdp;
        sce2c_c(LRO_CLOCK_ID, et, &sclkdp);

        // Try to get the spacecraft orientation at that time (spacecraft_from_J2000)
        SpiceDouble  spacecraft_from_J2000_R[3][3];
        SpiceDouble  clkout;
        SpiceBoolean orientationFound;
        ckgp_c(LRO_SPACECRAFT_CODE, sclkdp, tol, J2000_FRAME_STRING.c_str(), spacecraft_from_J2000_R, &clkout, &orientationFound);
        if (!orientationFound)
        {
          printf("SpiceEditor Error: Failed to obtain spacecraft position data at et %lf!!!!!!!!!!!!!!!!\n", et);
          return false;
        }

        // Try to get the planet orientation at that time (planet_from_J2000)
        SpiceDouble  planet_from_J2000_R[3][3];
        pxform_c(J2000_FRAME_STRING.c_str(), MOON_FRAME_STRING.c_str(), et, planet_from_J2000_R);
        // This matrix converts J2000 orientations to LRO orientations at et

        SpiceDouble  newPosition[3];
        if (localTransform) // Apply existing lronac position offset
        {
          // Convert the LRONAC offset from meters to kilometers
          SpiceDouble lronacOffsetKm[3];
          lronacOffsetKm[0] = lronacOffset[0] / 1000.0;
          lronacOffsetKm[1] = lronacOffset[1] / 1000.0;
          lronacOffsetKm[2] = lronacOffset[2] / 1000.0;
          
          SpiceDouble instOffset_J2000Frame[3];
          // Now need to convert the LRONAC offset (spacecraft to frame) into J2000 coordinates
          // - Multiply by the inverted rotation matrix to go from spacecraft frame to J2000 frame
          mtxv_c(spacecraft_from_J2000_R, lronacOffsetKm, instOffset_J2000Frame);
          for (int r=0; r<3; ++r) // Add the rotated offset to the original coordinate
            newPosition[r] = state[r] + instOffset_J2000Frame[r]; // Adding km to km here
        }
        else // Apply global transform from file
        {
          // Convert the planet transform to the correct coordinate space
          SpiceDouble fixedPlanet_from_J2000_R[3][3];
          SpiceDouble fixedPlanet_in_J2000_T[3];

          // Convert correcting rotation and translation into J2000 frame from planet frame
          mxm_c(planetFixed_from_planet_R, planet_from_J2000_R, fixedPlanet_from_J2000_R);

          // Convert the state into units of meters
          SpiceDouble stateMeters[3]; // J2000 frame
          for (int r=0; r<3; ++r)
            stateMeters[r] = state[r]*1000;

          // Apply the corrective planet-based transform to the vector currently represented in J2000 space
          SpiceDouble stateMoon[3], stateMoonFixed[3], stateJ2000Fixed[3];
          mxv_c(planet_from_J2000_R,       stateMeters, stateMoon);       // Convert to moon frame
          //mxv_c(planetFixed_from_planet_R, stateMoon,   stateMoonFixed);  // Apply correction --> Unused if transform is in local planet coordinates!
          for (int r=0; r<3; ++r) // Add in the translation
          {
//            double temp = stateMoonFixed[r];
            //stateMoonFixed[r] += planetFixed_from_planet_T[r]; 
            stateMoonFixed[r] = stateMoon[r] + planetFixed_from_planet_T[r]; // !AdjustedCameraModel applies transform about camera, not about the planet!
          }
          mtxv_c(planet_from_J2000_R, stateMoonFixed, stateJ2000Fixed); // Return to J2000 frame
          for (int r=0; r<3; ++r) // Convert from meters back to kilometers
          {
            newPosition[r] = stateJ2000Fixed[r] / 1000.0;
          }

        } // End of global transform case


        // Dump ET, state... line to a text file
        if (i > 0)
          outputFile << std::endl; // Add line breaks as needed
        outputFile.precision(16);
        outputFile << et;
        for (int q=0; q<3; ++q)
          outputFile << ", " << newPosition[q]; // Write out the new state
        for (int q=3; q<6; ++q)
          outputFile << ", " << state[q]; // Write input state from offset

      } // End of loop through steps

      outputFile.close(); // Close read data file


      //printf("Start time = %lf\n", startEt);
      //printf("Stop  time = %lf\n", lastEt);

    } // End loop through intervals for this body
    
    
  } // End of loop through bodies
  
  printf("Finished looping through kernel data, program finished.\n");

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

  return 0;
}






