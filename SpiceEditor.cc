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
#include "SpiceZfc.h"

#include <vw/Core/Exception.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Quaternion.h>

#include <asp/Tools/stereo.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageView.h>

//TODO: Would be nice to be able to remove these
using namespace vw;
using namespace vw::stereo;

using std::endl;
using std::setprecision;
using std::setw;

/* 
*/

struct Parameters : asp::BaseOptions 
{
  std::string inputPath;
  std::string outputPath;
  std::string leapSecondPath;
  //TODO: Add modification parameters
};


bool handle_arguments(int argc, char* argv[],
                     Parameters &opt) 
{ 
  po::options_description general_options("Options");
//  general_options.add_options()
//    ("intputPath", po::value(&opt.inputPath)->default_value(""), "Intput file path")
//    ("outputPath", po::value(&opt.outputPath)->default_value(""), "Output file path");
  
  general_options.add( asp::BaseOptionsDescription(opt) );
    
  po::options_description positional("");
  positional.add_options()
    ("inputPath",  po::value(&opt.inputPath))
    ("leapSecondPath",  po::value(&opt.leapSecondPath))
    ("outputPath", po::value(&opt.outputPath));  
    
  po::positional_options_description positional_desc;
  positional_desc.add("inputPath", 1);
  positional_desc.add("leapSecondPath", 1);
  positional_desc.add("outputPath", 1);
  

  std::string usage("[options] <Input Path> <Leap Second Path> <Output Path>");
  po::variables_map vm =
    asp::check_command_line( argc, argv, opt, general_options, general_options,
                             positional, positional_desc, usage );

  if ( !vm.count("inputPath") || !vm.count("leapSecondPath"))
    vw_throw( vw::ArgumentErr() << "Requires both inputs in order to proceed.\n\n"
              << usage << general_options );

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



//bool editSpiceFile(const std::string &inputPath, const std::string &leapSecondPath, const std::string &outputPath)
bool editSpiceFile(bool READ)  
{
  const std::string leapSecondPath;
  const std::string outputPath;

  /*
  //TODO: De-hardcode the extra input kernels
  // Load the input kernels
  furnsh_c(leapSecondPath.c_str());
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/base/kernels/pck/pck00009.tpc");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/base/kernels/spk/de118.bsp");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/base/kernels/spk/de245.bsp");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/base/kernels/spk/de405.bsp");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/base/kernels/lsk/kernels.0001.db");
  
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/pck/kernels.0002.db");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/tspk/kernels.0001.db");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/ck/kernels.1128.db");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/fk/kernels.0001.db");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/ik/kernels.0002.db");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/sclk/kernels.0004.db"); // Spacecraft clock
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/spk/kernels.1133.db");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/iak/kernels.0002.db");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/sclk/lro_clkcor_2013190_v00.tsc");
  
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/ik/lro_lroc_v16.ti");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/fk/lro_frames_2010277_v01.tf");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/fk/lro_frames_2012255_v02.tf");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/spk/fdf29_2010265_2010266_n01.bsp");
  */
  // Copied from ~/data/LinneCrater/M139829261LE.lronaccal.lronacecho.cub 
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/base/kernels/lsk/naif0010.tls");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/base/kernels/pck/pck00009.tpc");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/pck/moon_080317.tf");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/pck/moon_assoc_me.tf");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/tspk/moon_pa_de421_1900-2050.bpc");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/tspk/de421.bsp");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/ck/moc42_2010265_2010266_v01.bc");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/fk/lro_frames_2010277_v01.tf");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/ik/lro_lroc_v16.ti");
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/sclk/lro_clkcor_2013190_v00.tsc");
  
  std::string inputPath = "/home/smcmich1/testSpkFile.bsp";
  if (!READ)
    inputPath = "/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/spk/fdf29_2010265_2010266_n01.bsp";
  furnsh_c(inputPath.c_str());
  furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/iak/lro_instrumentAddendum_v03.ti");
  //furnsh_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/base/dems/ldem_128ppd_Mar2011_clon180_radius-_pad.cub");
  
  std::string frame   = "J2000";
  
//Reading kernelDbFile: /byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/base/dems/kernels.0005.db  
  
  //printf("Opened kernel %s\n", leapSecondPath.c_str());
  
  SpiceChar timeString[51];
  SPICEDOUBLE_CELL(cover, 2000);
  SpiceDouble b, e;  
  
  printf("CK file coverage:\n");
  // Check out the coverage of the ck (orientation) file
  SPICEINT_CELL(ckIds, 1000);
  ckobj_c ("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/ck/moc42_2010265_2010266_v01.bc", &ckIds);
  
  for (int i=0; i<card_c(&ckIds); ++i)
  {
    SpiceInt bodyId = SPICE_CELL_ELEM_I(&ckIds, i);
    
    // Get the coverage window
    SpiceBoolean  needav = 0;
    scard_c(0, &cover);
    ckcov_c("/byss/packages/isis-3.4.4-x86_64_linux_RHEL6/data/lro/kernels/ck/moc42_2010265_2010266_v01.bc",
            bodyId, needav, "SEGMENT", 0.0, "TDB", &cover);
       
   
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
    }
    
  }
  
  
    
//  furnsh_c(inputPath.c_str());  
//  printf("Opened kernel %s\n", inputPath.c_str());
 
  
  printf("SPK file coverage:\n");
  
  // Get list of bodies covered by the kernel
  SPICEINT_CELL(ids, 1000);
  spkobj_c(inputPath.c_str(), &ids);

  
  //TODO: Just loop through the times we pull from the tabledump
  
  // For each body
  for (int i=0; i<card_c(&ids); ++i)
  {
    SpiceInt bodyId = SPICE_CELL_ELEM_I(&ids, i);
    
    // Get the coverage window
    scard_c(0, &cover);
    spkcov_c(inputPath.c_str(), bodyId, &cover);

    
   
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

      // Set up data structures to hold modified data
      SpiceDouble *epochs          = new SpiceDouble[numSteps];
      SpiceDouble **modifiedStates = new SpiceDouble*[numSteps];
      for (int i=0; i<numSteps; ++i)
        modifiedStates[i] = new SpiceDouble[6];

      std::ofstream readFile;
      if (!READ)
        readFile.open("spkDataFile.txt");

      SpiceDouble lastEt;
      for (int i=0; i<numSteps; ++i)
      {

        SpiceDouble state[6];
        SpiceDouble lightTime;
        SpiceDouble et = startEt + stepSize*(SpiceDouble)i;
        epochs[i] = et;
        lastEt = et;
        
        //printf("Retrieving position at %lf\n", et);     
        
        std::string target  = "moon";  //TODO: Read these from file
        std::string observer= "earth";
        std::string absCorr = "NONE";
                
        // Retrieve the position of the spacecraft at this time
        //spkezr_c(target.c_str(),  et,  frame.c_str(), absCorr.c_str(), observer.c_str(), state, &lightTime);
        spkez_c(bodyId,  et,  frame.c_str(), absCorr.c_str(), 301, state, &lightTime); // Moon to LRO in J2000 frame

//        if (READ || (i % 10 == 0))
//          printf("ET %lf POS -->:  %lf, %lf, %lf -- %lf, %lf, %lf\n", et, state[0] * 1.0, state[1] * 1.0, state[2] * 1.0, state[3], state[4], state[5]); // Convert to meters
        
	
        SpiceInt     spacecraftClockId = -85;
        SpiceInt     inst = -85000; // Instrument ID code
        SpiceDouble  tol = 0; // Make sure we can do this 

        // Convert ephemeris time to spacecraft clock time
        SpiceDouble sclkdp;
        sce2c_c(spacecraftClockId, et, &sclkdp);

        // Try to get the spacecraft orientation at that time
        SpiceDouble  cmat[3][3];
        SpiceDouble  clkout;
        SpiceBoolean orientationFound;
        ckgp_c(-85000, sclkdp, tol, frame.c_str(), cmat, &clkout, &orientationFound);
        if (!orientationFound)
	      {
	        printf("Failed to obtain orientation data!!!!!!!!!!!!!!!!\n");
	        continue;
	      }
	
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

        // Now need to convert the LRONAC offsets (spacecraft to frame) into moon coordinates
        // - Multiply by the inverted rotation matrix to go from spacecraft frame to J2000 frame
        SpiceDouble  lronacOffsetL_J2000Frame[3];
        SpiceDouble  lronacOffsetR_J2000Frame[3];
        mtxvg_c(cmat, lronacOffsetL, 3, 3, lronacOffsetL_J2000Frame);
        mtxvg_c(cmat, lronacOffsetR, 3, 3, lronacOffsetR_J2000Frame);

        //printf("rotated offset = %lf, %lf, %lf\n", lronacOffsetL_J2000Frame[0], lronacOffsetL_J2000Frame[1], lronacOffsetL_J2000Frame[2]);
        //printf("\n");

        // Add the offsets to the camera position
        SpiceDouble*  modifiedState = modifiedStates[i];
        for (int n=0; n<3; ++n)
          modifiedState[n] = state[n] + lronacOffsetL_J2000Frame[n]; //TODO: Select correct offset!
//          modifiedState[n] = state[n]; // No change
        for (int n=3; n<6; ++n) // Copy the velocities with no change
          modifiedState[n] = state[n];


        if (!READ) // Dump ET, state... line to a text file
        {
          if (i > 0)
            readFile << std::endl; // Add line breaks as needed
          readFile.precision(16);
          readFile << et;
          for (int q=0; q<6; ++q)
            readFile << ", " << modifiedState[q];
        }

      } // End of loop through steps
      
      if (!READ)
        readFile.close(); // Close read data file


      printf("Start time = %lf\n", startEt);
      printf("Stop  time = %lf\n", lastEt);

      //if (!READ)
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
        */
        delete[] modifiedStates[i];
      }
      delete[] modifiedStates;

    } // End loop through intervals for this body
    
    
  } // End of loop through bodies
  
  printf("Finished looping through kernel data\n");
  
  // Unload the kernels
//  unload_c(leapSecondPath.c_str());
//  unload_c(inputPath.c_str());
  
  printf("Closed kernel\n");
  
  return true;
}




int main(int argc, char* argv[]) 
{
  try 
  {
    /*
    // Parse the input parameters
    Parameters params;
    if (!handle_arguments(argc, argv, params))
    {
      printf("Failed to parse input parameters!\n");
      return false;
    }
    */
    
    //editSpiceFile(params.inputPath, params.leapSecondPath, params.outputPath);
    
    printf("argc = %d\n", argc);
    bool read = (argc > 1);
    editSpiceFile(read);

  } ASP_STANDARD_CATCHES;

  return 0;
}






