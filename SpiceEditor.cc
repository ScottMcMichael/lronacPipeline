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

//#include <asp/Tools/stereo.h>
#include <stereo.h> // Using local version

//TODO: Would be nice to be able to remove these
using namespace vw;
//using namespace vw::stereo;

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
};


bool handle_arguments(int argc, char* argv[], Parameters &opt) 
{ 
  std::string outputPrefix;
  po::options_description general_options("Options");
  general_options.add_options()
    ("offsetCode",  po::value(&opt.offsetCode)->default_value(0), "Choose instrument: 0 = LE, 1 = RE")
    ("kernels", po::value<std::vector<std::string> >(&opt.kernelPaths)->multitoken(), "Paths to all required kernel files")
    ("transformFile", po::value<std::string>(&opt.transformFile)->default_value(""), "Path to pc_align global transform to apply")

    ("outputPrefix",  po::value(&outputPrefix)->default_value(""), "Output prefix");

  general_options.add( asp::BaseOptionsDescription(opt) );

  po::options_description positional("");
  positional.add_options();

  po::positional_options_description positional_desc;

  std::string usage("[options] ");
  po::variables_map vm =
    asp::check_command_line( argc, argv, opt, general_options, general_options,
                             positional, positional_desc, usage );

  if ( !vm.count("kernels") || !vm.count("outputPrefix"))
    vw_throw( vw::ArgumentErr() << "Requires kernels and output path in order to proceed.\n\n"
              << usage << general_options );

  // Finish paths to CK and SPK output data files
  opt.spkDataOutputPath = outputPrefix + "-spkData.txt";
  opt.ckDataOutputPath  = outputPrefix + "-ckData.txt";
  
  return true;
}

VW_DEFINE_EXCEPTION(SpiceErr, vw::Exception);
enum { LONG_MSG_LEN = 1840 };

void CHECK_SPICE_ERROR() 
{
  char longms[LONG_MSG_LEN + 1];

  longms[LONG_MSG_LEN] = 0;                // ensure the string is terminated
  if (failed_c()) 
  {
    getlms_(longms, LONG_MSG_LEN);
    reset_c();
    // Trim off the white space at the end
    for (int i = LONG_MSG_LEN-1; (i >= 0) && (longms[i] == ' '); --i)
      longms[i] = 0;
    std::cout
      << "SPICE: An error occured when accessing the SPICE information:\n\n"
      << longms;
    throw SpiceErr()
      << "SPICE: An error occured when accessing the SPICE information:\n\n"
      << longms;
  }
}


/// Loads a transform file from pc_align into a rotation and translation matrix
bool loadPcAlignTransform(const std::string &path, SpiceDouble R[][3], SpiceDouble *T)
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


//bool editSpiceFile(const std::string &inputPath, const std::string &leapSecondPath, const std::string &outputPath)
bool editSpiceFile(const Parameters &params)  
{
  const bool READ=false;

  const SpiceInt    MOON_CODE           = 301;
  const SpiceInt    LRO_CLOCK_ID        = -85;
  const SpiceInt    LRO_SPACECRAFT_CODE = -85000; // Instrument ID code
  const std::string J2000_FRAME_STRING  = "J2000";
  const std::string MOON_FRAME_STRING   = "IAU_MOON";

  //DEBUG
  // Set up georeference class with default moon datum
  vw::cartography::Datum datum("D_MOON");
  
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
  for (size_t i=0; i<params.kernelPaths.size(); ++i)
    furnsh_c(params.kernelPaths[i].c_str());
  
  SpiceDouble planetFixed_from_planet_R[3][3];
  SpiceDouble planetFixed_from_planet_T[3];
  if (!params.transformFile.empty())// Load the transform (it is in planet-space)
  {
    if (!loadPcAlignTransform(params.transformFile, planetFixed_from_planet_R, planetFixed_from_planet_T))
      return false;
  }
  
  SpiceChar timeString[51];
  SPICEDOUBLE_CELL(cover, 2000);
  SpiceDouble b, e;  
  
  // Check out the coverage of the ck (orientation) file
  printf("CK file coverage:\n");
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

    printf ( "\nCoverage for CK object %ld\n", (long int)bodyId );
    
    // Convert the coverage interval start and stop times to TDB calendar strings.
    for (int j=0; j<numIntervals; j++)
    {
      // Get the endpoints of the jth interval.
      wnfetd_c(&cover, j, &b, &e);

      // Convert the endpoints to TDB calendar format time strings and display them.
      timout_c(b, "YYYY MON DD HR:MN:SC.### (TDB) ::TDB", 51, timeString);

      printf( "\n"
              "Interval:  %ld\n"
	            "Start:     %s\n",
              j,
              timeString        );

      timout_c(e, "YYYY MON DD HR:MN:SC.### (TDB) ::TDB", 51, timeString);
      printf ( "Stop:      %s\n", timeString );
	
	    // Print raw times
      printf("b = %lf, e = %lf\n", b, e); 
    
      // If we are not applying a global transform there is nothing else to do here
      if (params.transformFile.empty())
        continue;
        
      // Correct the rotation at even intervals
      int numSteps = 5000;
      if (READ)
        numSteps = 20;
      SpiceDouble startEt = b+2; // Hack to avoid different start times
      SpiceDouble stepSize = (e - startEt) / numSteps;

    
      std::ofstream readFile;
      if (!READ)
      {
        printf("Writing file %s\n", params.ckDataOutputPath.c_str());
        readFile.open(params.ckDataOutputPath.c_str());
      }

      // For the number of specified steps
      SpiceDouble lastEt;
      for (int i=0; i<numSteps; ++i)
      {

        SpiceDouble state[6];
        SpiceDouble lightTime;
        SpiceDouble et = startEt + stepSize*(SpiceDouble)i;
        lastEt = et;
        
        std::string target   = "moon";  //TODO: Read these from file
        std::string observer = "earth";
        std::string absCorr  = "NONE";
        SpiceDouble  tol = 0; // Make sure we can do this 

        // Convert ephemeris time to spacecraft clock time
        SpiceDouble sclkdp;
        sce2c_c(LRO_CLOCK_ID, et, &sclkdp);

        //TODO: If possible move these transform requests to specific cases
        
        // Try to get the spacecraft orientation at that time (spacecraft_from_J2000)
        SpiceDouble  spacecraft_from_J2000_R[3][3];
        SpiceDouble  clkout;
        SpiceBoolean orientationFound;
        ckgp_c(LRO_SPACECRAFT_CODE, sclkdp, tol, J2000_FRAME_STRING.c_str(), spacecraft_from_J2000_R, &clkout, &orientationFound);
        if (!orientationFound)
        {
          printf("Failed to obtain spacecraft orientation data!!!!!!!!!!!!!!!!\n");
          return false;
        }
        /*
        printf("\n");
        printf("spacecraft_from_J2000_R:\n");
        for (int q=0; q<3; ++q)
        {
          for (int s=0; s<3; ++s)
          {
            printf("%lf, ", spacecraft_from_J2000_R[q][s]);
          }
          printf("\n");
        }
        */
        
        // Try to get the planet orientation at that time (planet_from_J2000)
        SpiceDouble  planet_from_J2000_R[3][3];
        pxform_c(J2000_FRAME_STRING.c_str(), MOON_FRAME_STRING.c_str(), et, planet_from_J2000_R);
        //if (!orientationFound)
        //{
        //  printf("Failed to obtain moon orientation data!!!!!!!!!!!!!!!!\n");
        //  return false;
        //}
        /*
        printf("planet_from_J2000_R:\n");
        for (int q=0; q<3; ++q)
        {
          for (int s=0; s<3; ++s)
          {
            printf("%lf, ", planet_from_J2000_R[q][s]);
          }
          printf("\n");
        }
        */
        
        // Convert the planet transform to the correct coordinate space
        SpiceDouble planetFixed_from_J2000_R[3][3];
     
        // Have fixedPlanet_from_planet (or reverse), sc and planet from J2000
        
        /*
        printf("planetFixed_from_planet_R:\n");
        for (int q=0; q<3; ++q)
        {
          for (int s=0; s<3; ++s)
          {
            printf("%lf, ", planetFixed_from_planet_R[q][s]);
          }
          printf("\n");
        }
        */
        
        // Convert correcting rotation and translation into J2000 frame from planet frame
        mxm_c(planetFixed_from_planet_R, planet_from_J2000_R, planetFixed_from_J2000_R);
        
        /*
        printf("planetFixed_from_J2000_R:\n");
        for (int q=0; q<3; ++q)
        {
          for (int s=0; s<3; ++s)
          {
            printf("%lf, ", planetFixed_from_J2000_R[q][s]);
          }
          printf("\n");
        }
        */
        
        // Apply the transform
        SpiceDouble spacecraftFixed_from_J2000_R[3][3];
        //mxm_c(planetFixed_from_J2000_R, spacecraft_from_J2000_R, spacecraftFixed_from_J2000_R); // Rotate the spacecraft orientation
       
        SpiceDouble spacecraft_from_Planet_R[3][3];
        SpiceDouble spacecraftFixed_from_Planet_R[3][3];
        mxmt_c(spacecraft_from_J2000_R,      planet_from_J2000_R,      spacecraft_from_Planet_R);       // Convert to moon frame
        mxm_c(planetFixed_from_planet_R,     spacecraft_from_Planet_R, spacecraftFixed_from_Planet_R);  // Apply correction
        mxm_c(spacecraftFixed_from_Planet_R, planet_from_J2000_R,      spacecraftFixed_from_J2000_R); // Return to J2000 frame
        
        
        /*
        printf("spacecraftFixed_from_J2000_R:\n");
        for (int q=0; q<3; ++q)
        {
          for (int s=0; s<3; ++s)
          {
            printf("%lf, ", spacecraftFixed_from_J2000_R[q][s]);
          }
          printf("\n");
        }
        */
        
        if (!READ) // Dump ET, new rotation... line to a text file
        {
          readFile.precision(16);
          readFile << et;
          for (int q=0; q<3; ++q)
            for (int s=0; s<3; ++s)
              readFile << " " << spacecraftFixed_from_J2000_R[q][s]; // Write out the new rotation
          readFile << std::endl;
        }
        
        
      } // End of steps loop
    
      if (!READ)
        readFile.close(); // Close read data file

    
    } // End of loop through intervals
    
  } // End of loop through bodies
  
  
  printf("SPK file coverage:\n");
  
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

    printf ( "\nCoverage for SPK object %ld\n", (long int)bodyId );
    
    // Convert the coverage interval start and stop times to TDB calendar strings.
    for (int j=0; j<numIntervals; j++)
    {
      // Get the endpoints of the jth interval.
      wnfetd_c(&cover, j, &b, &e);

      // Convert the endpoints to TDB calendar format time strings and display them.
      timout_c(b, "YYYY MON DD HR:MN:SC.### (TDB) ::TDB", 51, timeString);

      printf( "\n"
              "Interval:  %ld\n"
	            "Start:     %s\n",
              j,
              timeString        );

      timout_c(e, "YYYY MON DD HR:MN:SC.### (TDB) ::TDB", 51, timeString);
      printf ( "Stop:      %s\n", timeString );
	
	    // Print raw times
      printf("b = %lf, e = %lf\n", b, e); 
	
	
	    // Display the position at even intervals
      int numSteps = 5000;
      if (READ)
        numSteps = 20;
      SpiceDouble startEt = b+2; // Hack to avoid different start times
      SpiceDouble stepSize = (e - startEt) / numSteps;

      /*
      // Set up data structures to hold modified data
      SpiceDouble *epochs          = new SpiceDouble[numSteps];
      SpiceDouble **modifiedStates = new SpiceDouble*[numSteps];
      for (int i=0; i<numSteps; ++i)
        modifiedStates[i] = new SpiceDouble[6];
      */
  
      std::ofstream readFile;
      if (!READ)
      {
        printf("Writing file %s\n", params.spkDataOutputPath.c_str());
        readFile.open(params.spkDataOutputPath.c_str());
      }

      // For the number of specified steps
      SpiceDouble lastEt;
      for (int i=0; i<numSteps; ++i)
      {

        SpiceDouble state[6];
        SpiceDouble lightTime;
        SpiceDouble et = startEt + stepSize*(SpiceDouble)i;
        lastEt = et;
        
        //printf("Retrieving position at %lf\n", et);     
        
        std::string target  = "moon";  //TODO: Read these from file
        std::string observer= "earth";
        std::string absCorr = "NONE";
        const int MOON_CODE = 301;
        
        // Retrieve the position of the spacecraft at this time
        //spkezr_c(target.c_str(),  et,  frame.c_str(), absCorr.c_str(), observer.c_str(), state, &lightTime);
        //spkez_c(bodyId,  et,  frame.c_str(), absCorr.c_str(), 301, state, &lightTime); // Moon to LRO in J2000 frame
        spkez_c(bodyId,  et,  J2000_FRAME_STRING.c_str(), absCorr.c_str(), MOON_CODE, state, &lightTime); // LRO relative to Moon in J2000 frame, units are kilometers and km/sec

//        if (READ || (i % 10 == 0))
//          printf("ET %lf POS -->:  %lf, %lf, %lf -- %lf, %lf, %lf\n", et, state[0] * 1.0, state[1] * 1.0, state[2] * 1.0, state[3], state[4], state[5]); // Convert to meters
        
        SpiceDouble  tol = 0; // Make sure we can do this 

        // Convert ephemeris time to spacecraft clock time
        SpiceDouble sclkdp;
        sce2c_c(LRO_CLOCK_ID, et, &sclkdp);

        //TODO: If possible move these transform requests to specific cases
        
        // Try to get the spacecraft orientation at that time (spacecraft_from_J2000)
        SpiceDouble  spacecraft_from_J2000_R[3][3];
        SpiceDouble  clkout;
        SpiceBoolean orientationFound;
        ckgp_c(LRO_SPACECRAFT_CODE, sclkdp, tol, J2000_FRAME_STRING.c_str(), spacecraft_from_J2000_R, &clkout, &orientationFound);
        if (!orientationFound)
        {
          printf("Failed to obtain spacecraft orientation data!!!!!!!!!!!!!!!!\n");
          return false;
        }
        
        // Try to get the planet orientation at that time (planet_from_J2000)
        SpiceDouble  planet_from_J2000_R[3][3];
        pxform_c(J2000_FRAME_STRING.c_str(), MOON_FRAME_STRING.c_str(), et, planet_from_J2000_R);
        //if (!orientationFound)
        //{
        //  printf("Failed to obtain moon orientation data!!!!!!!!!!!!!!!!\n");
        //  return false;
        //}
        //TODO: Do we need the translation as well?
        
//        // TEST alternate method - seems to work exactly the same
//        SpiceDouble  planet_from_J2000_R_2[3][3];
//        tipbod_c(J2000_FRAME_STRING.c_str(), MOON_CODE, et, planet_from_J2000_R_2);
	
        // This matrix converts J2000 orientations to LRO orientations at et
/*        
        printf("ET %lf ROT -->:\n", et);
        printf("%lf, %lf, %lf\n", cmat[0][0], cmat[0][1], cmat[0][2]);
        printf("%lf, %lf, %lf\n", cmat[1][0], cmat[1][1], cmat[1][2]);
        printf("%lf, %lf, %lf\n", cmat[2][0], cmat[2][1], cmat[2][2]);
*/
/* I don't think we need this
        // Get rotation from spacecraft to camera frame
        SpiceDouble offsetRotation[3][3];
        pxform_c("LRO_SC_BUS", "LRO_LROCNACL", et, offsetRotation);
        printf("ET %lf INTERNAL ROT -->:\n", et);
        printf("%lf, %lf, %lf\n", offsetRotation[0][0], offsetRotation[0][1], offsetRotation[0][2]);
        printf("%lf, %lf, %lf\n", offsetRotation[1][0], offsetRotation[1][1], offsetRotation[1][2]);
        printf("%lf, %lf, %lf\n", offsetRotation[2][0], offsetRotation[2][1], offsetRotation[2][2]);
*/
        // Offsets from LRO spacecraft to the LRONAC cameras in kilometers
        SpiceDouble  lronacOffsetL[3] = {0.0013462, 0.0008890, -0.0001778}; 
        SpiceDouble  lronacOffsetR[3] = {0.0010160, 0.0008890, -0.0001778}; 

        
        SpiceDouble  newPosition[3];
        if (params.transformFile.empty()) // Apply existing lronac position offset
        {
          SpiceDouble  instOffset_J2000Frame[3];
          // Now need to convert the LRONAC offsets (spacecraft to frame) into J2000 coordinates
          // - Multiply by the inverted rotation matrix to go from spacecraft frame to J2000 frame
          if (params.offsetCode == 0)
            mtxv_c(spacecraft_from_J2000_R, lronacOffsetL, instOffset_J2000Frame);
          else
            mtxv_c(spacecraft_from_J2000_R, lronacOffsetR, instOffset_J2000Frame);
          for (int r=0; r<3; ++r) // Add the rotated offset to the original coordinate
            newPosition[r] = state[r] + instOffset_J2000Frame[r]; // Adding km to km here
        }
        else // Apply global transform from file
        {
          // Convert the planet transform to the correct coordinate space
          SpiceDouble fixedPlanet_from_J2000_R[3][3];
          SpiceDouble fixedPlanet_in_J2000_T[3];
  
          // Convert correcting rotation and translation into J2000 frame from planet frame
          mxm_c(planetFixed_from_planet_R, planet_from_J2000_R,       fixedPlanet_from_J2000_R);
          //mTxv_c(planet_from_J2000_R,      planetFixed_from_planet_T, fixedPlanet_in_J2000_T); 
          
          //printf("planetFixed_from_planet_T = %lf, %lf, %lf\n", planetFixed_from_planet_T[0], planetFixed_from_planet_T[1], planetFixed_from_planet_T[2]);
          //printf("fixedPlanet_in_J2000_T  = %lf, %lf, %lf\n", fixedPlanet_in_J2000_T[0], fixedPlanet_in_J2000_T[1], fixedPlanet_in_J2000_T[2]);
          /*
          printf("planetFixed_from_planet_R:\n");
          for (int q=0; q<3; ++q)
          {
            for (int s=0; s<3; ++s)
            {
              printf("%lf, ", planetFixed_from_planet_R[q][s]);
            }
            printf("\n");
          }
          printf("planet_from_J2000_R:\n");
          for (int q=0; q<3; ++q)
          {
            for (int s=0; s<3; ++s)
            {
              printf("%lf, ", planet_from_J2000_R[q][s]);
            }
            printf("\n");
          }
          printf("fixedPlanet_from_J2000_R:\n");
          for (int q=0; q<3; ++q)
          {
            for (int s=0; s<3; ++s)
            {
              printf("%lf, ", fixedPlanet_from_J2000_R[q][s]);
            }
            printf("\n");
          }
          */
          
          // Convert the state into units of meters
          SpiceDouble stateMeters[3]; // J2000 frame
          for (int r=0; r<3; ++r)
            stateMeters[r] = state[r]*1000;
          
          // Apply the corrective planet-based transform to the vector currently represented in J2000 space
          SpiceDouble stateMoon[3], stateMoonFixed[3], stateJ2000Fixed[3];
          mxv_c(planet_from_J2000_R,       stateMeters, stateMoon);       // Convert to moon frame
          mxv_c(planetFixed_from_planet_R, stateMoon,   stateMoonFixed);  // Apply correction
          for (int r=0; r<3; ++r) // Add in the translation
          {
            double temp = stateMoonFixed[r];
            stateMoonFixed[r] += planetFixed_from_planet_T[r];
            //printf("stateMoon[%d]: %lf --> %lf + %lf = %lf\n", r, stateMoon[r], temp, planetFixed_from_planet_T[r], stateMoonFixed[r]);
          }
          mtxv_c(planet_from_J2000_R, stateMoonFixed, stateJ2000Fixed); // Return to J2000 frame
          for (int r=0; r<3; ++r) // Convert from meters back to kilometers
          {
            newPosition[r] = stateJ2000Fixed[r] / 1000.0;
            //printf("stateMeters[%d]: %lf --> %lf\n", r, stateMeters[r], newPosition[r]*1000.0);
          }
/*
          mxv_c(fixedPlanet_from_J2000_R, stateMeters, newPosition); // Rotate the point about the origin
          for (int r=0; r<3; ++r) // Add in the translation
          {
            //newPosition[r] = state[r]; //DEBUG
            double temp = newPosition[r];
            newPosition[r] += fixedPlanet_in_J2000_T[r];
            printf("stateMeters[%d]: %lf --> %lf + %lf = %lf\n", r, stateMeters[r], temp, fixedPlanet_in_J2000_T[r], newPosition[r]);
            newPosition[r] = newPosition[r] / 1000.0; // Convert the positions back into kilometers
          }
*/
/*
          // Convert from GCC to GDC
          Vector3 gccCoordOld(stateMoon[0],      stateMoon[1],      stateMoon[2]);
          Vector3 gccCoordNew(stateMoonFixed[0], stateMoonFixed[1], stateMoonFixed[2]);
          Vector3 gdcCoordOld = datum.cartesian_to_geodetic(gccCoordOld);
          Vector3 gdcCoordNew = datum.cartesian_to_geodetic(gccCoordNew);
          std::cout << "Old GDC = " << gdcCoordOld << std::endl;
          std::cout << "New GDC = " << gdcCoordNew << std::endl;
          printf("\n");
*/
          
        } // End of global transform case
          
        //printf("rotated offset = %lf, %lf, %lf\n", lronacOffsetL_J2000Frame[0], lronacOffsetL_J2000Frame[1], lronacOffsetL_J2000Frame[2]);
        //printf("\n");

	/*
        // Add the offsets to the camera position
        SpiceDouble*  modifiedState = modifiedStates[i];
        for (int n=0; n<3; ++n)
          modifiedState[n] = state[n] + lronacOffsetL_J2000Frame[n]; //TODO: Select correct offset!
//          modifiedState[n] = state[n]; // No change
        for (int n=3; n<6; ++n) // Copy the velocities with no change
          modifiedState[n] = state[n];
	*/

        if (!READ) // Dump ET, state... line to a text file
        {
          if (i > 0)
            readFile << std::endl; // Add line breaks as needed
          readFile.precision(16);
          readFile << et;
          for (int q=0; q<3; ++q)
            readFile << ", " << newPosition[q]; // Write out the new state
          for (int q=3; q<6; ++q)
            readFile << ", " << state[q]; // Write input state from offset
        }

      } // End of loop through steps
      
      if (!READ)
        readFile.close(); // Close read data file


      printf("Start time = %lf\n", startEt);
      printf("Stop  time = %lf\n", lastEt);

      /*
      if (false)
      {
        // Create a new SPK file
        SpiceInt    newSpkFileHandle; 
        std::string spkfileOutputPath = "/home/smcmich1/testSpkFile.bsp";
        SpiceInt    reservedCommentSize=2000;
        spkopn_c(spkfileOutputPath.c_str(), "GENERATED_SPK_FILE", reservedCommentSize, &newSpkFileHandle);


        // Now save the adjusted position to an output CK file
        // SPK type 13, polynomial degree 11
        SpiceInt degree = 11; // This must be odd
        SpiceInt center = 301; // MOON
        SpiceInt nSteps = numSteps;
        spkw13_c(newSpkFileHandle, bodyId, 301, frame.c_str(), startEt, lastEt,
                    "Segment_Identifier", // Anything to put here?
                    degree, nSteps, modifiedStates, epochs);
                    
        // Close the new SPK file
        spkcls_c(newSpkFileHandle);
      }
      */  

      /*
      printf("---------\n");
      // Delete the allocated data
      delete[] epochs;
      for (int i=0; i<numSteps; ++i)
      {
        /*
        if (!READ && (i % 50 == 0))
        {
          for (int m=0; m<6; ++m)
            printf("%lf, ", modifiedStates[i][m]); // print written values
          delete[] modifiedStates[i];
          printf("\n");
        }
        else
        delete[] modifiedStates[i];
      }
      delete[] modifiedStates;
      */

    } // End loop through intervals for this body
    
    
  } // End of loop through bodies
  
  printf("Finished looping through kernel data\n");
  
  // Unload the kernels
  for (size_t i=0; i<params.kernelPaths.size(); ++i)
    unload_c(params.kernelPaths[i].c_str());  
    
  printf("Closed kernel\n");
  
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
    
    /*
    // DEBUG: Set fields manually
    Parameters params;
    params.offsetCode = 0;
    params.outputPath = "/home/smcmich1"
    params.kernelPaths.push_back(("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/base/kernels/lsk/naif0010.tls");
    params.kernelPaths.push_back(("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/base/kernels/pck/pck00009.tpc");
    params.kernelPaths.push_back(("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/pck/moon_080317.tf");
    params.kernelPaths.push_back(("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/pck/moon_assoc_me.tf");
    params.kernelPaths.push_back(("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/tspk/moon_pa_de421_1900-2050.bpc");
    params.kernelPaths.push_back(("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/tspk/de421.bsp");
    params.kernelPaths.push_back(("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/ck/moc42_2010265_2010266_v01.bc");
    params.kernelPaths.push_back(("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/fk/lro_frames_2010277_v01.tf");
    params.kernelPaths.push_back(("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/ik/lro_lroc_v16.ti");
    params.kernelPaths.push_back(("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/sclk/lro_clkcor_2013190_v00.tsc");
    //std::string inputPath = "/home/smcmich1/testSpkFile.bsp";
    params.kernelPaths.push_back("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/spk/fdf29_2010265_2010266_n01.bsp");
    params.kernelPaths.push_back("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/iak/lro_instrumentAddendum_v03.ti");
    //furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/base/dems/ldem_128ppd_Mar2011_clon180_radius-_pad.cub");    
    */

    editSpiceFile(params);

  } ASP_STANDARD_CATCHES;

  return 0;
}






