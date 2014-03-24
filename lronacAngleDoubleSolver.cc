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


/// \file lronacAngleDoubleSolver.cc
///

#include <iostream>
#include <iomanip>

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
#include <asp/IsisIO/IsisInterfaceLineScan.h>

#include <stereo.h>

#include <lronacSolverSupport.h>
#include <lronacSolverModelDouble.h>

using namespace vw;
using namespace vw::stereo;
using namespace asp::isis;
using namespace asp;
using std::endl;
using std::setprecision;
using std::setw;



struct Parameters : asp::BaseOptions 
{
  // Input paths
  std::string leftFilePath;
  std::string rightFilePath;
  std::string leftStereoFilePath;
  std::string rightStereoFilePath;
  
  
  std::string outputPrefix;
  
  std::string matchingLeftPointsPath; 
  std::string matchingRightPointsPath;
  std::string matchingLeftCrossPointsPath;
  std::string matchingRightCrossPointsPath;
  
  bool initialOnly; ///< If true only compute starting state, don't run the solver.
  
  std::string initialValuePath;

  double expectedSurfaceElevation;
  int cropWidth; ///< Specifies image overlap for use with ipfind
  
  bool debug;
};



bool handle_arguments(int argc, char* argv[],
                     Parameters &opt) 
{ 
  po::options_description general_options("Options");
  general_options.add_options()
    ("debug",                  po::bool_switch(&opt.debug                 )->default_value(false),  "DEBUG mode")
    ("outputPrefix",                 po::value      (&opt.outputPrefix                )->required(),            "Output prefix to use")
    ("initialOnly",                  po::bool_switch(&opt.initialOnly                 )->default_value(false),  "Just compute initial state (don't solve)")
    ("initialValues",                po::value      (&opt.initialValuePath            )->default_value(""),     "Path to file containing state parameter values (probably from previous output)")
    ("crop-width",                   po::value      (&opt.cropWidth                   )->default_value(200),    "Crop images to this width before disparity search")
    ("elevation",                    po::value      (&opt.expectedSurfaceElevation    )->default_value(0.0),    "Start solver estimate at this surface elevation")
    ("matchingPixelsLeftPath",       po::value      (&opt.matchingLeftPointsPath      )->default_value(""),     "Path to left-leftS stereo pixel file")
    ("matchingPixelsRightPath",      po::value      (&opt.matchingRightPointsPath     )->default_value(""),     "Path to right-rightS stereo pixel file")
    ("matchingPixelsLeftCrossPath",  po::value      (&opt.matchingLeftCrossPointsPath )->default_value(""),     "Path to left-rightS stereo pixel file")
    ("matchingPixelsRightCrossPath", po::value      (&opt.matchingRightCrossPointsPath)->default_value(""),     "Path to right-leftS stereo pixel file")
    ("leftCubePath",                 po::value      (&opt.leftFilePath                )->default_value(""),     "Path to left input cube")
    ("rightCubePath",                po::value      (&opt.rightFilePath               )->default_value(""),     "Path to right input cube")
    ("leftStereoCubePath",           po::value      (&opt.leftStereoFilePath          )->default_value(""),     "Path to left stereo input cube")
    ("rightStereoCubePath",          po::value      (&opt.rightStereoFilePath         )->default_value(""),     "Path to right stereo input cube");
  
  general_options.add( asp::BaseOptionsDescription(opt) );
    
  po::options_description positional("");
    
  po::positional_options_description positional_desc;

  std::string usage("[options]");
  po::variables_map vm =
    asp::check_command_line( argc, argv, opt, general_options, general_options,
                             positional, positional_desc, usage );

//  if ()
//    vw_throw( ArgumentErr() << "Missing required input arguments!.\n\n"
//              << usage << general_options );

  return true;
}

//-------------------------------------------------------------------------------------------

// Load matching points from a file
bool loadMatchingPixels(const std::string &pointPath, PointObsList &pixelVals)
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
    //printf("p: %d, --> %lf, %lf, %lf, %lf\n", p, sample1, line1, sample2, line2);
    ++p;
  }
  file.close();
  
  // Convert from std::vector to packed pixel objects
  pixelVals.leftObsList.resize (p);
  pixelVals.rightObsList.resize(p);
  for (int i=0; i<p; ++i)
  {
    pixelVals.leftObsList [i]  = Vector2(col1[i], row1[i]);
    pixelVals.rightObsList[i] = Vector2(col2[i], row2[i]);
  }
  
  return true;
}

//-------------------------------------------------------------------------------------------


// Search for matching pixels in the LE/RE overlap
bool findMatchingPixels(const std::string &leftFilePath, const std::string &rightFilePath, const std::string &logFilePath, const int overlapWidth, PointObsList &pixelVals)
{
 
  // Load both images  
  printf("Loading images left=%s and right=%s...\n",
         leftFilePath.c_str(),
         rightFilePath.c_str());
  DiskImageView<PixelGray<float> > left_disk_image (leftFilePath );
  DiskImageView<PixelGray<float> > right_disk_image(rightFilePath);
  
  printf("Left  input image size: %d rows, %d cols\n", left_disk_image.rows(),  left_disk_image.cols());
  printf("Right input image size: %d rows, %d cols\n", right_disk_image.rows(), right_disk_image.cols());
  
  //const int imageWidth      = std::min(left_disk_image.cols(), right_disk_image.cols());
  const int imageHeight     = std::min(left_disk_image.rows(), right_disk_image.rows());
  const int imageTopRow     = 0;
  //const int imageMidPointX  = imageWidth / 2;
  //const int cropStartX      = imageMidPointX - (params.cropWidth/2);

  // Restrict processing to the border of the images
  // - Since both images were nproj'd the overlap areas should be in about the same spots.
  const int leftStartX = left_disk_image.cols()-overlapWidth;
  const BBox2i leftRoi (leftStartX, imageTopRow, overlapWidth, imageHeight);
  const BBox2i rightRoi(0,          imageTopRow, overlapWidth, imageHeight);
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
  double inlierThreshold     = 10.0; // Want to be somewhat generous here
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

  // This value is governed by how slow the solver gets with large numbers of points
  // TODO: Fiddle with values to get this number higher?
  const int TARGET_NUM_POINTS = 200;

  // Convert the matching points into the correct format
  const int pointSkip = ransac_ip1.size() / TARGET_NUM_POINTS; // Pick skip to get about the desired number of points

  // Initialize random seed to generate random starting offset
  srand(time(0));
  const int startIndex = rand() % pointSkip;
  printf("Starting index = %d, point skip = %d\n", startIndex, pointSkip);
  
  const size_t numMatchedPts = floor(ransac_ip1.size() / (pointSkip+1));
  printf("Num sampled points = %lu\n", numMatchedPts);
  //leftRow.set_size(numMatchedPts), leftCol.set_size(numMatchedPts), rightRow.set_size(numMatchedPts), rightCol.set_size(numMatchedPts);
  //int i = startIndex;
  /*
  std::ofstream pointFile("/home/smcmich1/initialPixelDiff.csv");
  for (size_t p=0; p<numMatchedPts; ++p)
  {
    leftCol [p] = ransac_ip1[i][0] + leftStartX;
    leftRow [p] = ransac_ip1[i][1];
    rightCol[p] = ransac_ip2[i][0];
    rightRow[p] = ransac_ip2[i][1];
    printf("p: %d, i: %d --> %lf, %lf, %lf, %lf\n", p, i, leftCol[p], leftRow[p], rightCol[p], rightRow[p]);
    pointFile << leftRow[p] - rightRow[p] << ", " << ransac_ip1[i][0] - rightCol[p] << std::endl;
    i+= pointSkip; 
  }
  pointFile.close();
  */
  
  // Convert from std::vector to packed pixel objects
  // - Also record matching pixel files
  printf("Recording IpFind matches to %s\n", logFilePath.c_str());
  std::ofstream outputFile(logFilePath.c_str());
  pixelVals.leftObsList.resize (numMatchedPts);
  pixelVals.rightObsList.resize(numMatchedPts);
  int i = startIndex;
  for (int p=0; p<numMatchedPts; ++p)
  {
    pixelVals.leftObsList [p] = Vector2(ransac_ip1[i][0] + leftStartX, ransac_ip1[i][1]);
    pixelVals.rightObsList[p] = Vector2(ransac_ip2[i][0],              ransac_ip2[i][1]);
    i+= pointSkip; 
    outputFile << pixelVals.leftObsList[p][0] << ", " << pixelVals.leftObsList[p][1] << ", " << pixelVals.rightObsList[p][0] << ", " << pixelVals.rightObsList[p][1] << std::endl;
  }
  outputFile.close();
  
  return true;
}

//-------------------------------------------------------------------------------------------

/// Load all available points based on the input parameters
/// - Returns the number of points loaded or zero if there is an error
size_t loadInputPointPairs(const Parameters &params, PointObsList &leftPixelPairs,      PointObsList &rightPixelPairs,
                                                     PointObsList &overlapPairs,        PointObsList &stereoOverlapPairs,
                                                     PointObsList &leftCrossPixelPairs, PointObsList &rightCrossPixelPairs)
{

  // If these files already exist read them in, otherwise compute them and save them.
  std::string mainIpFindPath   = params.outputPrefix + "-mainIpFindPixels.csv";
  std::string stereoIpFindPath = params.outputPrefix + "-stereoIpFindPixels.csv";

  // -- Adjacent cubes section (ipfind based matches) --
  if (params.initialOnly == false) // Don't do any of this if we are only computing initial state
  {
    if (boost::filesystem::exists(boost::filesystem::path(mainIpFindPath)))
    {
      // Point file already exists
      printf("Loading list of matched pixels from file %s\n", mainIpFindPath.c_str());
      if (!loadMatchingPixels(mainIpFindPath, overlapPairs))
        return 0;
    }
    else // Point file does not exist
    {
      // If the left and right cubes were passed in, call ipfind in the overlap region to get points.
      if ((params.leftFilePath.size() > 0) && (params.rightFilePath.size() > 0))
      {
        printf("Searching for matching pixels in image overlap region\n");
        if (!findMatchingPixels(params.leftFilePath, params.rightFilePath, mainIpFindPath, params.cropWidth, overlapPairs))
          return 0;
      }
    }
    if (boost::filesystem::exists(boost::filesystem::path(stereoIpFindPath)))
    {
      // Point file already exists
      printf("Loading list of matched pixels from file %s\n", stereoIpFindPath.c_str());
      if (!loadMatchingPixels(stereoIpFindPath, stereoOverlapPairs))
        return 0;
    }
    else // Point file does not exist
    {
      // If the left and right stereo cubes were passed in, call ipfind in the overlap region to get points.
      if ((params.leftStereoFilePath.size() > 0) && (params.rightStereoFilePath.size() > 0))
      {
        printf("Searching for matching pixels in stereo image overlap region\n");
        if (!findMatchingPixels(params.leftStereoFilePath, params.rightStereoFilePath, stereoIpFindPath, params.cropWidth, stereoOverlapPairs))
          return 0;
      }
    }
  } // End adjacent cubes section

  // -- Stereo cubes section (stereo based matches) --
  if (params.matchingLeftPointsPath.size() > 0)
  {
    printf("Loading list of matched pixels from file %s\n", params.matchingLeftPointsPath.c_str());
    if (!loadMatchingPixels(params.matchingLeftPointsPath, leftPixelPairs))
      return 0;
  }
  if (params.matchingRightPointsPath.size() > 0)
  {
    printf("Loading list of matched pixels from file %s\n", params.matchingRightPointsPath.c_str());
    if (!loadMatchingPixels(params.matchingRightPointsPath, rightPixelPairs))
      return 0;
  }
  if (params.matchingLeftCrossPointsPath.size() > 0)
  {
    printf("Loading list of matched pixels from file %s\n", params.matchingLeftCrossPointsPath.c_str());
    if (!loadMatchingPixels(params.matchingLeftCrossPointsPath, leftCrossPixelPairs))
      return false;
  }
  if (params.matchingRightCrossPointsPath.size() > 0)
  {
    printf("Loading list of matched pixels from file %s\n", params.matchingRightCrossPointsPath.c_str());
    if (!loadMatchingPixels(params.matchingRightCrossPointsPath, rightCrossPixelPairs))
      return false;
  }

  // Total up all the loaded points
  const size_t totalNumPoints = overlapPairs.size()        + stereoOverlapPairs.size() +
                                leftPixelPairs.size()      + rightPixelPairs.size() +
                                leftCrossPixelPairs.size() + rightCrossPixelPairs.size();
  return totalNumPoints;
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

// Main solver function
bool optimizeRotations(Parameters & params)
{
  // This is the number of non-point parameters there are.
  int NUM_CAMERA_PARAMS = 12;

  printf("**************************************************************************\n");
  printf("*** Started LRONAC double angle solver.\n");
  
  // Get boost paths to images, used to check if they exist
  boost::filesystem::path leftBoostPath       (params.leftFilePath );
  boost::filesystem::path rightBoostPath      (params.rightFilePath);
  boost::filesystem::path leftStereoBoostPath (params.leftStereoFilePath);
  boost::filesystem::path rightStereoBoostPath(params.rightStereoFilePath);
  
  // Load each of the four input cubes into the solver if they are present
  LrocPairModel lrocClass;
  if (params.leftFilePath.size() > 0)
  {
    if (boost::filesystem::exists(boost::filesystem::path(params.leftFilePath)))
      lrocClass.loadLeftCamera(params.leftFilePath);
    else
    {
      printf("Error: input file %s is missing!\n", params.leftFilePath.c_str());
      return false;
    }
  }
  else // No file given, make sure we don't need it!
  {
    if ( (params.matchingLeftPointsPath.size()      > 0) ||
         (params.matchingLeftCrossPointsPath.size() > 0)  )
    {
      printf("Error: Specified input file requires left input image!\n");
      return false;
    }
  }

  if (params.rightFilePath.size() > 0)
  {
    if (boost::filesystem::exists(boost::filesystem::path(params.rightFilePath)))
      lrocClass.loadRightCamera(params.rightFilePath);
    else
    {
      printf("Error: input file %s is missing!\n", params.rightFilePath.c_str());
      return false;
    }
  }
  else // No file given, make sure we don't need it!
  {
    if ( (params.matchingRightPointsPath.size()      > 0) ||
         (params.matchingRightCrossPointsPath.size() > 0)  )
    {
      printf("Error: Specified input file requires right input image!\n");
      return false;
    }
  }

  if (params.leftStereoFilePath.size() > 0)
  {
    if (boost::filesystem::exists(boost::filesystem::path(params.leftStereoFilePath)))
      lrocClass.loadLeftStereoCamera(params.leftStereoFilePath);
    else
    {
      printf("Error: input file %s is missing!\n", params.leftStereoFilePath.c_str());
      return false;
    }
  }
  else // No file given, make sure we don't need it!
  {
    if ( (params.matchingLeftPointsPath.size()       > 0) ||
         (params.matchingRightCrossPointsPath.size() > 0)  )
    {
      printf("Error: Specified input file requires left stereo input image!\n");
      return false;
    }
  }

  if (params.rightStereoFilePath.size() > 0)
  {
    if (boost::filesystem::exists(boost::filesystem::path(params.rightStereoFilePath)))
      lrocClass.loadRightStereoCamera(params.rightStereoFilePath);
    else
    {
      printf("Error: input file %s is missing!\n", params.rightStereoFilePath.c_str());
      return false;
    }
  }
  else // No file given, make sure we don't need it!
  {
    if ( (params.matchingRightPointsPath.size()     > 0) ||
         (params.matchingLeftCrossPointsPath.size() > 0)  )
    {
      printf("Error: Specified input file requires right stereo input image!\n");
      return false;
    }
  }
  // Finished loading cubes into the camera model


  if ((params.rightFilePath.size() > 0) && (!boost::filesystem::exists(boost::filesystem::path(params.rightFilePath))))
  {
    printf("Error: input file %s is missing!\n", params.rightFilePath.c_str());
    return false;
  }
  if ((params.leftStereoFilePath.size() > 0) && (!boost::filesystem::exists(boost::filesystem::path(params.leftStereoFilePath))))
  {
    printf("Error: input file %s is missing!\n", params.leftStereoFilePath.c_str());
    return false;
  }
  if ((params.rightStereoFilePath.size() > 0) && (!boost::filesystem::exists(boost::filesystem::path(params.rightStereoFilePath))))
  {
    printf("Error: input file %s is missing!\n", params.rightStereoFilePath.c_str());
    return false;
  }

 
  // Now load up all of the pixel pairs
  PointObsList leftPixelPairs, rightPixelPairs, overlapPairs, stereoOverlapPairs, leftCrossPixelPairs, rightCrossPixelPairs;
  const size_t totalNumPoints = loadInputPointPairs(params, leftPixelPairs,      rightPixelPairs,
                                                            overlapPairs,        stereoOverlapPairs,
                                                            leftCrossPixelPairs, rightCrossPixelPairs);
  if (totalNumPoints <= 40)
  {
    printf("Error: Did not load enough points to compute a solution!");
    return false;
  }
  
  //TODO: Can't really load point values without a rotation-modified pixel to vector function.
  // If a path to an initial value file was provided, load them
  size_t estimatedNumParams = NUM_CAMERA_PARAMS;// + totalNumPoints*PARAMS_PER_POINT;
  std::vector<double> initialValues;
  if (!params.initialValuePath.empty())
  {
    initialValues.reserve(estimatedNumParams);
    printf("Reading initial values from file %s\n", params.initialValuePath.c_str());
    std::ifstream initalValueFile(params.initialValuePath.c_str());
    size_t initialValueCount = 0;
    //while (initialValueFile.good())
    for (int i=0; i<estimatedNumParams; ++i) 
    {
      double newVal;
      initalValueFile >> newVal;
      initialValues.push_back(newVal);
      ++initialValueCount;
    }
    initalValueFile.close();
    if ((initalValueFile.fail()) || (initialValueCount > estimatedNumParams))
    {
      printf("Error reading from initial value file %s\n", params.initialValuePath.c_str());
      printf("Read %d values, expected %d values\n", initialValueCount, estimatedNumParams);
      return false;
    }
  }

  
  
  // Load the inital points into the solver
  printf("Initializing solver state...\n");
  Vector<double> initialState;
  std::vector<double> initialErrorMeters;
  if (!lrocClass.getInitialStateEstimate(overlapPairs,        stereoOverlapPairs,
                                         leftPixelPairs,      rightPixelPairs,
                                         leftCrossPixelPairs, rightCrossPixelPairs,
                                         initialState,        initialErrorMeters,
                                         params.expectedSurfaceElevation, initialValues))
  {
    printf("Failed to get initial state!\n");
    return false;
  }
  // The initial state contains the camera parameters, then the point coordinates of each set of points in sequence

  // Write initial state error to file
  std::string initialErrorMetersPath = params.outputPrefix + "-initialPointErrorMeters.csv";
  std::ofstream initialErrorMetersFile(initialErrorMetersPath.c_str());
  for (size_t i=0; i<initialErrorMeters.size(); ++i)
    initialErrorMetersFile << initialErrorMeters[i] << std::endl;
  initialErrorMetersFile.close();


  // Set up georeference class with default moon datum
  vw::cartography::Datum datum("D_MOON");

  ////printf("Writing initial state log to %s\n", initialStatePath.c_str());
  //std::ofstream initialObsFile("/home/smcmich1/initialObsVector.csv");
  //for (size_t i=0; i<packedObservations.size(); ++i)
  //  initialObsFile << packedObservations[i] << std::endl;
  //initialObsFile.close();
  
  //return false;
  
  //Vector<double> initialPredictions = lrocClass(initialState);
  //std::ofstream initialPredFile("/home/smcmich1/initialObsPrediction.csv");
  //for (size_t i=0; i<initialPredictions.size(); ++i)
  //  initialPredFile << initialPredictions[i] << std::endl;
  //initialPredFile.close();

  std::string initialStatePath = params.outputPrefix + "-initialParamState.csv";
  printf("Writing initial state log to %s\n", initialStatePath.c_str());
  std::ofstream initialStateFile(initialStatePath.c_str());
  for (size_t i=0; i<initialState.size(); ++i)
    initialStateFile << initialState[i] << std::endl;
  initialStateFile.close();


  // Write initial points as GDC coordinates for google earth
  std::string initialGdcCoordPath = params.outputPrefix + "-initialGdcPoints.csv";
  std::ofstream initialGdcCoordFile;
  initialGdcCoordFile.open(initialGdcCoordPath.c_str());
  for (size_t i=0; i<totalNumPoints; ++i)
  {
    // Convert from GCC to GDC
    Vector3 gccPoint(initialState[NUM_CAMERA_PARAMS+ 3*i], initialState[NUM_CAMERA_PARAMS+ 3*i+1], initialState[NUM_CAMERA_PARAMS+ 3*i+2]);
    Vector3 gdcCoord = datum.cartesian_to_geodetic(gccPoint);
    initialGdcCoordFile.precision(12);
    // Write lat, lon, height so pc_align tool can read these files
    initialGdcCoordFile << gdcCoord[1] << ", " << gdcCoord[0] << ", " << gdcCoord[2] << std::endl;
  }
  initialGdcCoordFile.close();


  // Compute the initial error - euclidean point distance
  std::string initialErrorPath = params.outputPrefix + "-initialPointError.csv";
  printf("Writing initial error log to %s\n", initialErrorPath.c_str());
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
  printf(">>>> Mean point error before optimization = %lf <<<<<<\n", meanInitialError);


  //Vector<double> initialComputedObservations = lrocClass(initialState);
  
  if (params.initialOnly)
  {
    printf("Stopping after initial state calculation\n");
    return true;
  }
  
  //@@@@@@@@@@@@@@@@@@@@ SOLVER BLOCK @@@@@@@@@@@@@@@@@@@@@@@@@
  printf("Running solver...\n");
  Vector<double> finalParams;
  printf("Using Ceres solver!\n");
  //TODO: Turn on Glog
  
  // Create Ceres solver object
  ceres::Problem problem;
  
  // Set up camera parameters for solver
  double* localRotation       = &(initialState[0]);
  double* globalRotation      = &(initialState[3]);
  double* globalPosition      = &(initialState[6]);
  double* localStereoRotation = &(initialState[9]);

  problem.AddParameterBlock(localRotation,       3);
  problem.AddParameterBlock(globalRotation,      3);
  problem.AddParameterBlock(globalPosition,      3);
  problem.AddParameterBlock(localStereoRotation, 3);

  const int NUM_PARAMS_CAMERA_SECTION  = 3;
  const int NUM_PARAMS_PER_POINT       = 3;
  const int NUM_PARAMS_PER_OBSERVATION = 2;

  size_t currentPointIndex = NUM_CAMERA_PARAMS;
  ceres::LossFunction* lossFunction = new ceres::CauchyLoss(5.0);

  if (overlapPairs.size() > 0)
    printf("Loading parameters for main camera pair...\n");
  for (size_t i=0; i<overlapPairs.size(); ++i) // For each input point
  {
    // Set up this point's parameters for the solver
    double* pointParams = &(initialState[currentPointIndex]);
    //printf("Loading L-R point %lf, %lf, %lf\n", pointParams[0], pointParams[1], pointParams[2]);
    currentPointIndex += NUM_PARAMS_PER_POINT;
    problem.AddParameterBlock(pointParams, NUM_PARAMS_PER_POINT);
    
    // Add the function and residual block for the left camera
    ceres::CostFunction* costFunctionLeft = 
            new ceres::NumericDiffCostFunction<LeftCostFunctor, ceres::CENTRAL, NUM_PARAMS_PER_OBSERVATION, NUM_PARAMS_PER_POINT>(
                  new LeftCostFunctor(&lrocClass, overlapPairs.leftObsList[i]));

    problem.AddResidualBlock(costFunctionLeft, lossFunction, pointParams);

    // Add the function and residual block for the right camera (slightly more complex)
    ceres::CostFunction* costFunctionRight = 
            new ceres::NumericDiffCostFunction<RightCostFunctor, ceres::CENTRAL, NUM_PARAMS_PER_OBSERVATION, NUM_PARAMS_CAMERA_SECTION, NUM_PARAMS_PER_POINT>(
                  new RightCostFunctor(&lrocClass, overlapPairs.rightObsList[i]));
    
    problem.AddResidualBlock(costFunctionRight, lossFunction, localRotation, pointParams);
    
  } // End of loop through main camera pair points


  //currentPointIndex += overlapPairs.size()*NUM_PARAMS_PER_POINT; //DEBUG!!!!!!

  if (stereoOverlapPairs.size() > 0)
    printf("Loading parameters for stereo camera pair...\n");
  for (size_t i=0; i<stereoOverlapPairs.size(); ++i) // For each input point
  {
    // Set up this point's parameters for the solver
    double* pointParams = &(initialState[currentPointIndex]);
    currentPointIndex += NUM_PARAMS_PER_POINT;
    problem.AddParameterBlock(pointParams, NUM_PARAMS_PER_POINT);
  
    
    // Add the function and residual block for the left camera
    ceres::CostFunction* costFunctionLeft = 
            new ceres::NumericDiffCostFunction<LeftStereoCostFunctor, ceres::CENTRAL, NUM_PARAMS_PER_OBSERVATION, NUM_PARAMS_CAMERA_SECTION, NUM_PARAMS_CAMERA_SECTION, NUM_PARAMS_PER_POINT>(
                  new LeftStereoCostFunctor(&lrocClass, stereoOverlapPairs.leftObsList[i]));

    problem.AddResidualBlock(costFunctionLeft, lossFunction, globalRotation, globalPosition, pointParams);

    // Add the function and residual block for the right camera (slightly more complex)
    ceres::CostFunction* costFunctionRight = 
            new ceres::NumericDiffCostFunction<RightStereoCostFunctor, ceres::CENTRAL, NUM_PARAMS_PER_OBSERVATION, NUM_PARAMS_CAMERA_SECTION, NUM_PARAMS_CAMERA_SECTION, NUM_PARAMS_CAMERA_SECTION, NUM_PARAMS_PER_POINT>(
                  new RightStereoCostFunctor(&lrocClass, stereoOverlapPairs.rightObsList[i]));
    
    problem.AddResidualBlock(costFunctionRight, lossFunction, globalRotation, globalPosition, localStereoRotation, pointParams);
    
  } // End of loop through stereo camera pair points


  //currentPointIndex += stereoOverlapPairs.size()*NUM_PARAMS_PER_POINT; //DEBUG!!!!!!

  if (leftPixelPairs.size() > 0)
    printf("Loading parameters for two left cameras...\n");
  for (size_t i=0; i<leftPixelPairs.size(); ++i) // For each input point
  {
    // Set up this point's parameters for the solver
    double* pointParams = &(initialState[currentPointIndex]);
    //printf("Loading L-L point %lf, %lf, %lf\n", pointParams[0], pointParams[1], pointParams[2]);
    currentPointIndex += NUM_PARAMS_PER_POINT;
    problem.AddParameterBlock(pointParams, NUM_PARAMS_PER_POINT);
    
    // Add the function and residual block for the left camera
    ceres::CostFunction* costFunctionLeft = 
            new ceres::NumericDiffCostFunction<LeftCostFunctor, ceres::CENTRAL, NUM_PARAMS_PER_OBSERVATION, NUM_PARAMS_PER_POINT>(
                  new LeftCostFunctor(&lrocClass, leftPixelPairs.leftObsList[i]));

    problem.AddResidualBlock(costFunctionLeft, lossFunction, pointParams);


    // Add the function and residual block for the stereo left camera (slightly more complex)
    ceres::CostFunction* costFunctionRight = 
            new ceres::NumericDiffCostFunction<LeftStereoCostFunctor, ceres::CENTRAL, NUM_PARAMS_PER_OBSERVATION, NUM_PARAMS_CAMERA_SECTION, NUM_PARAMS_CAMERA_SECTION, NUM_PARAMS_PER_POINT>(
                  new LeftStereoCostFunctor(&lrocClass, leftPixelPairs.rightObsList[i]));
    
    problem.AddResidualBlock(costFunctionRight, lossFunction, globalRotation, globalPosition, pointParams);
    
  } // End of loop through both left camera points


  //currentPointIndex += leftPixelPairs.size()*NUM_PARAMS_PER_POINT; //DEBUG!!!!!!

  if (rightPixelPairs.size() > 0)
    printf("Loading parameters for two right cameras...\n");
  for (size_t i=0; i<rightPixelPairs.size(); ++i) // For each input point
  {
    // Set up this point's parameters for the solver
    double* pointParams = &(initialState[currentPointIndex]);
    currentPointIndex += NUM_PARAMS_PER_POINT;
    problem.AddParameterBlock(pointParams, NUM_PARAMS_PER_POINT);
  
    
    // Add the function and residual block for the right camera
    ceres::CostFunction* costFunctionLeft = 
            new ceres::NumericDiffCostFunction<RightCostFunctor, ceres::CENTRAL, NUM_PARAMS_PER_OBSERVATION, NUM_PARAMS_CAMERA_SECTION, NUM_PARAMS_PER_POINT>(
                  new RightCostFunctor(&lrocClass, rightPixelPairs.leftObsList[i]));

    problem.AddResidualBlock(costFunctionLeft, lossFunction, localRotation, pointParams);

    // Add the function and residual block for the stereo right camera (slightly more complex)
    ceres::CostFunction* costFunctionRight = 
            new ceres::NumericDiffCostFunction<RightStereoCostFunctor, ceres::CENTRAL, NUM_PARAMS_PER_OBSERVATION, NUM_PARAMS_CAMERA_SECTION, NUM_PARAMS_CAMERA_SECTION, NUM_PARAMS_CAMERA_SECTION, NUM_PARAMS_PER_POINT>(
                  new RightStereoCostFunctor(&lrocClass, rightPixelPairs.rightObsList[i]));
    
    problem.AddResidualBlock(costFunctionRight, lossFunction, globalRotation, globalPosition, localStereoRotation, pointParams);

  } // End of loop through both right camera points


  if (leftCrossPixelPairs.size() > 0)
    printf("Loading parameters for left cross camera pair...\n");
  for (size_t i=0; i<leftCrossPixelPairs.size(); ++i) // For each input point
  {
    // Set up this point's parameters for the solver
    double* pointParams = &(initialState[currentPointIndex]);
    //printf("Loading L-R point %lf, %lf, %lf\n", pointParams[0], pointParams[1], pointParams[2]);
    currentPointIndex += NUM_PARAMS_PER_POINT;
    problem.AddParameterBlock(pointParams, NUM_PARAMS_PER_POINT);

    // Add the function and residual block for the left camera
    ceres::CostFunction* costFunctionLeft =
            new ceres::NumericDiffCostFunction<LeftCostFunctor, ceres::CENTRAL, NUM_PARAMS_PER_OBSERVATION, NUM_PARAMS_PER_POINT>(
                  new LeftCostFunctor(&lrocClass, leftCrossPixelPairs.leftObsList[i]));

    problem.AddResidualBlock(costFunctionLeft, lossFunction, pointParams);

    // Add the function and residual block for the stereo right camera (slightly more complex)
    ceres::CostFunction* costFunctionRight =
            new ceres::NumericDiffCostFunction<RightStereoCostFunctor, ceres::CENTRAL, NUM_PARAMS_PER_OBSERVATION, NUM_PARAMS_CAMERA_SECTION, NUM_PARAMS_CAMERA_SECTION, NUM_PARAMS_CAMERA_SECTION, NUM_PARAMS_PER_POINT>(
                  new RightStereoCostFunctor(&lrocClass, leftCrossPixelPairs.rightObsList[i]));

    problem.AddResidualBlock(costFunctionRight, lossFunction, globalRotation, globalPosition, localStereoRotation, pointParams);

  } // End of loop through left cross camera pair points

  if (rightCrossPixelPairs.size() > 0)
    printf("Loading parameters for right cross camera pair...\n");
  for (size_t i=0; i<rightCrossPixelPairs.size(); ++i) // For each input point
  {
    // Set up this point's parameters for the solver
    double* pointParams = &(initialState[currentPointIndex]);
    //printf("Loading L-R point %lf, %lf, %lf\n", pointParams[0], pointParams[1], pointParams[2]);
    currentPointIndex += NUM_PARAMS_PER_POINT;
    problem.AddParameterBlock(pointParams, NUM_PARAMS_PER_POINT);

    // Add the function and residual block for the left camera
    ceres::CostFunction* costFunctionLeft =
            new ceres::NumericDiffCostFunction<LeftStereoCostFunctor, ceres::CENTRAL, NUM_PARAMS_PER_OBSERVATION, NUM_PARAMS_CAMERA_SECTION, NUM_PARAMS_CAMERA_SECTION, NUM_PARAMS_PER_POINT>(
                  new LeftStereoCostFunctor(&lrocClass, rightCrossPixelPairs.leftObsList[i]));

    problem.AddResidualBlock(costFunctionLeft, lossFunction, globalRotation, globalPosition, pointParams);

    // Add the function and residual block for the right camera
    ceres::CostFunction* costFunctionRight =
            new ceres::NumericDiffCostFunction<RightCostFunctor, ceres::CENTRAL, NUM_PARAMS_PER_OBSERVATION, NUM_PARAMS_CAMERA_SECTION, NUM_PARAMS_PER_POINT>(
                  new RightCostFunctor(&lrocClass, rightCrossPixelPairs.rightObsList[i]));

    problem.AddResidualBlock(costFunctionRight, lossFunction, localRotation, pointParams);

  } // End of loop through right cross camera pair points


  printf("Finished loading points into the solver!\n");
       
  // TODO: Select solver options
  ceres::Solver::Options solverOptions;
  solverOptions.max_num_iterations           = 50; //TODO: Play with these again!
  //solverOptions.linear_solver_type           = ceres::DENSE_SCHUR;
  solverOptions.linear_solver_type           = ceres::SPARSE_SCHUR;
  solverOptions.minimizer_progress_to_stdout = true;
  solverOptions.max_num_line_search_direction_restarts = 8;
  solverOptions.use_nonmonotonic_steps = false; // Allow non-descent steps to try to find global minimum --> Seems to lead to bad results!
  solverOptions.max_num_consecutive_invalid_steps = 10;
  solverOptions.num_threads = 1; // SPICE cannot handle multiple threads here!
  solverOptions.num_linear_solver_threads = 8; 
  solverOptions.use_approximate_eigenvalue_bfgs_scaling = false; // Use approximate eigenvalue scaling
  //solverOptions.solver_log = "~/data/ceresOutput.txt";
  // There are many more options to play with!
  
  // Execute the Ceres solver
  printf("Starting the Ceres solver...\n");
  ceres::Solver::Summary summary;
  ceres::Solve(solverOptions, &problem, &summary);
  
  //std::ofstream ceresLog("/home/smcmich1/data/ceresOutput2.txt");
  std::cout << summary.FullReport() << "\n";
  //ceresLog << summary.FullReport();
  //ceresLog.close();

  printf("All done with the Ceres solver!\n");

  printf("Copying out Ceres solver results.\n");
  // Get the information back from the solver!
  finalParams.set_size(initialState.size());
  for (size_t i=0; i<initialState.size(); ++i)
    finalParams[i] = initialState[i];
  
  //return true;
  

  //@@@@@@@@@@@@@@@@@@@@@@@@ END SOLVER BLOCK @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  std::string finalStatePath = params.outputPrefix + "-finalParamState.csv";
  printf("Writing final state log to %s\n", finalStatePath.c_str());
  std::ofstream finalStateFile(finalStatePath.c_str());
  for (size_t i=0; i<finalParams.size(); ++i)
    finalStateFile << finalParams[i] << std::endl;
  finalStateFile.close();


  // The computeError() function calculates the mean euclidean distance from the observation points of the LE and RE cameras for each point.
  std::string finalErrorPath = params.outputPrefix + "-finalPointError.csv";
  printf("Writing final error log to %s\n", finalErrorPath.c_str());

  std::ofstream finalErrorFile(finalErrorPath.c_str());
  double meanFinalError = 0;
  std::vector<double> finalError = lrocClass.computeError(finalParams);
  //Vector<double> finalPredictions = lrocClass(finalParams);
  //Vector<double> rawError(finalPredictions.size());
  //std::ofstream predictionFile("/home/smcmich1/finalPredictions.csv");
  //for (size_t i=0; i<finalPredictions.size(); ++i)
  //{
  //  predictionFile << finalPredictions[i] << std::endl;
    //rawError[i] = packedObservations[i] - finalPredictions[i];
    //finalErrorFile << packedObservations[i] - finalPredictions[i] << std::endl;
  //  if (i % 4 == 0)
  //    meanFinalError += finalError[i/4];
  //}

  for (size_t i=0; i<finalError.size(); ++i)
  {
    finalErrorFile << finalError[i] << std::endl;
    meanFinalError += finalError[i];
  }
  meanFinalError = meanFinalError / finalError.size();

  //predictionFile.close();
  finalErrorFile.close();
  //meanFinalError = meanFinalError / currentError.size();
  
  // Compute the median error
  std::sort(finalError.begin(), finalError.end());
  double medianError = 0;
  size_t centralIndex = finalError.size()/2;
  if ((finalError.size() % 2) == 0)
    medianError = (finalError[centralIndex-1] + finalError[centralIndex]) / 2.0;
  else
    medianError = finalError[centralIndex];
  printf(">>>> Median final error = %lf <<<<<<\n", medianError);


  // Write initial points as GDC coordinates for google earth
  std::string outputGdcPath = params.outputPrefix + "-outputGdcPoints.csv";
  std::ofstream finalGdcCoordFile(outputGdcPath.c_str());
  for (size_t i=0; i<totalNumPoints; ++i)
  {
    // Convert from GCC to GDC
    Vector3 gccPoint(finalParams[NUM_CAMERA_PARAMS+ 3*i], finalParams[NUM_CAMERA_PARAMS+ 3*i+1], finalParams[NUM_CAMERA_PARAMS+ 3*i+2]);
    Vector3 gdcCoord = datum.cartesian_to_geodetic(gccPoint);
    finalGdcCoordFile.precision(12);
    // Write lat, lon, height so pc_align tool can read these files
    finalGdcCoordFile << gdcCoord[1] << ", " << gdcCoord[0] << ", " << gdcCoord[2] << std::endl;
  }
  finalGdcCoordFile.close();

  // Compare the intersection locations of the points to the solved locations of the points
  vw::Vector<double> intersectedParams;
  std::vector<double> finalIntersectErrorMeters, finalParamsCopy;
  finalParamsCopy.resize(initialState.size());
  for (size_t i=0; i<initialState.size(); ++i)
    finalParamsCopy[i] = finalParams[i];
  lrocClass.getInitialStateEstimate(overlapPairs,        stereoOverlapPairs,
                                    leftPixelPairs,      rightPixelPairs,
                                    leftCrossPixelPairs, rightCrossPixelPairs,
                                    intersectedParams,   finalIntersectErrorMeters,
                                    params.expectedSurfaceElevation, finalParamsCopy);

  std::string gdcComparePath = params.outputPrefix + "-outputGdcPointsIntersected.csv";
  std::ofstream outputGdcIntersected( gdcComparePath.c_str() );
  for (size_t i=0; i<totalNumPoints; ++i)
  {
    // Convert from GCC to GDC
    Vector3 gccPoint(intersectedParams[NUM_CAMERA_PARAMS+ 3*i  ],
                     intersectedParams[NUM_CAMERA_PARAMS+ 3*i+1],
                     intersectedParams[NUM_CAMERA_PARAMS+ 3*i+2]);
    Vector3 gdcCoord = datum.cartesian_to_geodetic(gccPoint);
    outputGdcIntersected.precision(12);
    // Write lat, lon, height so pc_align tool can read these files
    outputGdcIntersected << gdcCoord[1] << ", " << gdcCoord[0] << ", " << gdcCoord[2] << std::endl;
  }
  outputGdcIntersected.close();




  printf(">>>> Mean point error after optimization = %lf <<<<<<\n", meanFinalError);
  printf(">>>> Mean error change = %lf <<<<<<\n", meanFinalError - meanInitialError);
  
  // --------------- Summary of results ------------------------------
  const double rad2deg = 180.0 / M_PI;
  const double deg2rad = M_PI / 180.0;
  printf("Local  rotation angles (degrees): %lf, %lf, %lf\n", finalParams[0]*rad2deg, finalParams[1]*rad2deg, finalParams[2]*rad2deg);
  printf("Local  rotation angles (radians): %lf, %lf, %lf\n", finalParams[0], finalParams[1], finalParams[2]);
  printf("Global rotation angles (degrees): %lf, %lf, %lf\n", finalParams[3]*rad2deg, finalParams[4]*rad2deg, finalParams[5]*rad2deg);
  printf("Global rotation angles (radians): %lf, %lf, %lf\n", finalParams[3], finalParams[4], finalParams[5]);
  printf("Global translation     (meters ): %lf, %lf, %lf\n", finalParams[6], finalParams[7], finalParams[8]);
  printf("Stereo local rotation  (degrees): %lf, %lf, %lf\n", finalParams[9]*rad2deg, finalParams[10]*rad2deg, finalParams[11]*rad2deg);
  printf("Stereo local rotation  (radians): %lf, %lf, %lf\n", finalParams[9], finalParams[10], finalParams[11]);
  
  // ----------------- Now handle the local rotation ------------------------
  vw::math::Matrix<double,3,3> frameRotation = vw::math::euler_to_rotation_matrix(1.129*deg2rad, -0.079*deg2rad, 180.0*deg2rad, "xyz");
  std:: cout << "Base LRONAC inter-camera rotation angles (degrees) = Rx(1.129) * Ry(-0.079) * Rz(180.0)"  << std::endl;

  // Existing offset matrix (from the fk file) = sc_from_instrument
  
  // This is the rotation matrix we solved for (matches code in IsisInterfaceLineScan.cc)
  // solvedRotation = instrument(good)_from_instrument
  vw::math::Matrix<double,3,3> tempRot = vw::math::euler_to_rotation_matrix(finalParams[0], finalParams[1], finalParams[2], "xyz");
  vw::math::Matrix<double,3,3> localRotationMatrix = transpose(frameRotation * transpose(tempRot)); // == sc_from_instrument(good)
  std::cout << "Solved local rotation matrix: " << localRotationMatrix << std::endl;

  // Extract rotation angles --> Returns Rx(a)*Ry(b)*Rz(c)
  // -- Angle order needs to be  1, 2, 3
  Vector3 outputLocalAngles = vw::math::rotation_matrix_to_euler_xyz(localRotationMatrix);
  // Output angles (rZ, rY, rX)
  printf("Combined local rotation angles (degrees): %lf, %lf, %lf\n", outputLocalAngles[2]*rad2deg, outputLocalAngles[1]*rad2deg, outputLocalAngles[0]*rad2deg);
  
  // --------------- Handle global rotation (translation is easy) -----------------------
  vw::math::Matrix<double,3,3> globalRotMatrix = vw::math::euler_to_rotation_matrix(finalParams[3], finalParams[4], finalParams[5], "xyz");
  std::cout << "Global rotation matrix: " << globalRotMatrix << std::endl;


  // --------------- Now handle the stereo local rotation ----------------------
  tempRot = vw::math::euler_to_rotation_matrix(finalParams[9], finalParams[10], finalParams[11], "xyz");
  vw::math::Matrix<double,3,3> stereoLocalRotationMatrix = transpose(frameRotation * transpose(tempRot)); // == sc_from_instrument(good)
  std::cout << "Solved stereo local rotation matrix: " << stereoLocalRotationMatrix << std::endl;

  // Extract rotation angles --> Returns Rx(a)*Ry(b)*Rz(c)
  // -- Angle order needs to be  1, 2, 3
  Vector3 outputStereoLocalAngles = vw::math::rotation_matrix_to_euler_xyz(stereoLocalRotationMatrix);
  // Output angles (rZ, rY, rX)
  printf("Combined stereo local rotation angles (degrees): %lf, %lf, %lf\n", outputStereoLocalAngles[2]*rad2deg, outputStereoLocalAngles[1]*rad2deg, outputStereoLocalAngles[0]*rad2deg);


  // Three output files are generated,
  std::string localRotationPath       = params.outputPrefix + "-localRotationMatrix.csv";
  std::string globalTransformPath     = params.outputPrefix + "-globalTransformMatrix.csv";
  std::string stereoLocalRotationPath = params.outputPrefix + "-stereoLocalRotationMatrix.csv";

  // Write out the local rotation
  printf("Writing output file %s\n", localRotationPath.c_str());
  std::ofstream outputFile(localRotationPath.c_str());
  for (unsigned int r=0; r<localRotationMatrix.rows(); ++r)
  {
    for (unsigned int c=0; c<localRotationMatrix.cols(); ++c)
    {
      outputFile << localRotationMatrix(r,c) << std::endl;
    }
  }
  outputFile.close();


  // Write out the global rotation/translation matrix that was applied to the camera
  printf("Writing output file %s\n", globalTransformPath.c_str());
  outputFile.open(globalTransformPath.c_str());
  vw::math::Vector<double,3>   axisAngle   (finalParams[3], finalParams[4], finalParams[5]);
  vw::math::Quaternion<double> rotQuat     (axis_angle_to_quaternion(axisAngle));
  vw::math::Matrix<double,3,3> globalRotMat(rotQuat.rotation_matrix());

  for (int q=0; q<3; ++q)
    outputFile << globalRotMat[q][0] << " " << globalRotMat[q][1] << " " << globalRotMat[q][2] << " " << finalParams[6+q] << std::endl;
  outputFile << "0 0 0 1" << std::endl;
  outputFile.close();

  // Write out the stereo local rotation
  printf("Writing output file %s\n", stereoLocalRotationPath.c_str());
  outputFile.open(stereoLocalRotationPath.c_str());
  for (unsigned int r=0; r<stereoLocalRotationMatrix.rows(); ++r)
  {
    for (unsigned int c=0; c<stereoLocalRotationMatrix.cols(); ++c)
    {
      outputFile << stereoLocalRotationMatrix(r,c) << std::endl;
    }
  }
  outputFile.close();

  printf("**** LRONAC double angle solver finished!\n");
  printf("**************************************************************************\n");

  return true;
}


//-------------------------------------------------------------------------------------------

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

    if (params.debug)
    {
      IsisInterfaceLineScanRot leftStereoCam(params.leftStereoFilePath);
      AdjustedCameraModelRot rotCam(boost::shared_ptr<IsisInterfaceLineScanRot>(&leftStereoCam, boost::serialization::null_deleter()));
      Vector3 rotVec(0.00323453, -0.0112794, -0.0145248);
      Vector3 offsetVec(3.22023, -4.80386, 23.8804);
      rotCam.set_axis_angle_rotation(rotVec);
      rotCam.set_translation(offsetVec);

      Vector3 nullVec(0,0,0);
      int line   = 500;
      int sample = 500;
      Vector2 pixel(sample, line);
      double et = 318930609.91952;
      printf("et = %lf\n", et);

      Matrix3x3 instMatrix, bodyMatrix;
      leftStereoCam.getMatricesAtTime(et, instMatrix, bodyMatrix);
      for (int r=0; r<3; ++r)
      {
        for (int c=0; c<3; ++c)
        {
          printf("%lf, ", instMatrix[r][c]);
        }
        printf("\n");
      }

      // Input pixel pair: 3200,800    3430,997
      // Computed GDC coordinate: 31.3108323222, -67.1141237031, 255.905626238
      // Computed GCC coordinate: 577350, -1.36772e+06, 903026


      /*
      Vector3 point(577350, -1.36772e+06, 903026);
      Vector2 pixel, idealPixel(3430,997);
      pixel = rotCam.point_to_pixel(point);
      std::cout << "pixel = " << pixel << std::endl;
      std::cout << "pixel diff = " << pixel  - idealPixel << std::endl;

      pixel = rotCam.point_to_pixel_rotated(point, nullVec, 997);
      std::cout << "pixel 2= " << pixel << std::endl;
      std::cout << "pixel 2 diff = " << pixel  - idealPixel << std::endl;
      */
      return 0;
    }

    
    optimizeRotations(params);
  } ASP_STANDARD_CATCHES;

  return 0;
}



