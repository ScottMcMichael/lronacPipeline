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


/// \file lronacAngleSolver.cc
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

#include <lronacSolverSupport.h>
#include <lronacSolverModel.h>

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
  std::string outputPath;
  std::string gdcPointsOutPath;
  std::string matchingPointsPath; ///< If not set, ipfind/ransac is used to find matching points

  bool worldTransform;  ///< Solve for transform in world coordinate frame instead of camera frame
  bool includePosition; ///< Also solve for a position offset in the specified coordinate frame

  bool initialOnly; ///< If true only compute starting state, don't run the solver.
  
  //std::vector<double> initialValues; ///< Starting state parameters
  std::string initialValuePath;

  int  cropWidth; ///< Specifies image overlap for use with ipfind
};



bool handle_arguments(int argc, char* argv[],
                     Parameters &opt) 
{ 
  po::options_description general_options("Options");
  general_options.add_options()
    ("outputPath",         po::value      (&opt.outputPath        )->default_value(""),     "Write angles to this path (in degrees)")
    ("matchingPixelsPath", po::value      (&opt.matchingPointsPath)->default_value(""),     "File to load matching points from")
    ("gdcPointsOutPath",   po::value      (&opt.gdcPointsOutPath  )->default_value(""),     "Output path for final GDC points")
    ("worldTransform",     po::bool_switch(&opt.worldTransform    )->default_value(false),  "Compute transform in world frame instead of camera frame")
    ("includePosition",    po::bool_switch(&opt.includePosition   )->default_value(false),  "Also solve for a 3D position offset (in chosen frame)")
    ("initialOnly",        po::bool_switch(&opt.initialOnly       )->default_value(false),  "Just compute initial state (don't solve)")
    ("initialValues",      po::value      (&opt.initialValuePath  )->default_value(""),     "Path to file containing state parameter values (probably from previous output)")
    ("crop-width",         po::value      (&opt.cropWidth         )->default_value(200),    "Crop images to this width before disparity search");
  
  general_options.add( asp::BaseOptionsDescription(opt) );
    
  po::options_description positional("");
  positional.add_options()
    ("left",  po::value(&opt.leftFilePath))
    ("right", po::value(&opt.rightFilePath));  
    
  po::positional_options_description positional_desc;
  positional_desc.add("left",  1);
  positional_desc.add("right", 1);
  

  std::string usage("[options] <left> <right>");
  po::variables_map vm =
    asp::check_command_line( argc, argv, opt, general_options, general_options,
                             positional, positional_desc, usage );

  if ( !vm.count("left") || !vm.count("right") )
    vw_throw( ArgumentErr() << "Requires <left> and <right> input in order to proceed.\n\n"
              << usage << general_options );

  return true;
}

//-------------------------------------------------------------------------------------------

// Load mathing points from a file
bool loadMatchingPixels(const std::string &pointPath,
          Vector<double> &leftRow, Vector<double> &leftCol, Vector<double> &rightRow, Vector<double> &rightCol)
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
  
  leftRow.set_size (p);
  leftCol.set_size (p);
  rightRow.set_size(p);
  rightCol.set_size(p);
  for (int i=0; i<p; ++i)
  {
    leftCol [i] = col1[i];
    leftRow [i] = row1[i];
    rightCol[i] = col2[i];
    rightRow[i] = row2[i];
  }
  
  return true;
}

//-------------------------------------------------------------------------------------------

// Search for matching pixels in the LE/RE overlap
bool findMatchingPixels(const Parameters &params, Vector<double> &leftRow, Vector<double> &leftCol, Vector<double> &rightRow, Vector<double> &rightCol)
{
 
  // Load both images  
  printf("Loading images left=%s and right=%s...\n",
         params.leftFilePath.c_str(),
         params.rightFilePath.c_str());
  DiskImageView<PixelGray<float> > left_disk_image (params.leftFilePath );
  DiskImageView<PixelGray<float> > right_disk_image(params.rightFilePath);
  
  printf("Left  input image size: %d rows, %d cols\n", left_disk_image.rows(),  left_disk_image.cols());
  printf("Right input image size: %d rows, %d cols\n", right_disk_image.rows(), right_disk_image.cols());
  
  //const int imageWidth      = std::min(left_disk_image.cols(), right_disk_image.cols());
  const int imageHeight     = std::min(left_disk_image.rows(), right_disk_image.rows());
  const int imageTopRow     = 0;
  //const int imageMidPointX  = imageWidth / 2;
  //const int cropStartX      = imageMidPointX - (params.cropWidth/2);

  // Restrict processing to the border of the images
  // - Since both images were nproj'd the overlap areas should be in about the same spots.
  const int leftStartX = left_disk_image.cols()-params.cropWidth;
  const BBox2i leftRoi (leftStartX, imageTopRow, params.cropWidth, imageHeight);
  const BBox2i rightRoi(0,          imageTopRow, params.cropWidth, imageHeight);
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
  double inlierThreshold     = 5.0; // Want to be somewhat generous here
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
  const int TARGET_NUM_POINTS = 70;

  // Initialize random seed, then generate secret number between 0 and TARGET_NUM_POINTS-1
  srand(time(0));
  const int startIndex = rand() % TARGET_NUM_POINTS;
  printf("Starting index = %d\n", startIndex);

  // Convert the matching points into the correct format
  const int    pointSkip     = ransac_ip1.size() / TARGET_NUM_POINTS; // Pick skip to get the right number of points
  const size_t numMatchedPts = ransac_ip1.size() / pointSkip;
  printf("Num sampled points = %lu\n", numMatchedPts);
  leftRow.set_size(numMatchedPts), leftCol.set_size(numMatchedPts), rightRow.set_size(numMatchedPts), rightCol.set_size(numMatchedPts);
  int i = startIndex;
  for (size_t p=0; p<numMatchedPts; ++p)
  {
    leftCol [p] = ransac_ip1[i][0] + leftStartX;
    leftRow [p] = ransac_ip1[i][1];
    rightCol[p] = ransac_ip2[i][0];
    rightRow[p] = ransac_ip2[i][1];
    //printf("p: %d, i: %d --> %lf, %lf, %lf, %lf\n", p, i, leftCol[p], leftRow[p], rightCol[p], rightRow[p]);
    i+= pointSkip; 
  }  
  
  return true;
}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

// Main solver function
bool optimizeRotations(Parameters & params)
{
  
  // Verify images are present
  boost::filesystem::path leftBoostPath (params.leftFilePath );
  boost::filesystem::path rightBoostPath(params.rightFilePath);
  
  if (!boost::filesystem::exists(boost::filesystem::path(params.leftFilePath)))
  {
    printf("Error: input file %s is missing!\n", params.leftFilePath.c_str());
    return false;
  }
  if (!boost::filesystem::exists(boost::filesystem::path(params.rightFilePath)))
  {
    printf("Error: input file %s is missing!\n", params.rightFilePath.c_str());
    return false;
  }

  // How many non-point parameters are there?
  int startOfPts = 3;
  if (params.includePosition)
    startOfPts = 6;

  // If a path to an initial value file was provided, load them
  std::vector<double> initialValues;
  if (!params.initialValuePath.empty())
  {
    initialValues.resize(startOfPts);
    printf("Reading inital values from file %s\n", params.initialValuePath.c_str());
    std::ifstream initalValueFile(params.initialValuePath.c_str());
    for (int i=0; i<startOfPts; ++i)
    {
      double newVal;
      initalValueFile >> newVal;
      initialValues[i] = newVal;
    }
    initalValueFile.close();
    if (initalValueFile.fail())
    {
      printf("Error reading from inital value file %s\n", params.initialValuePath.c_str());
      return false;
    }
  }

  printf("Constructing geometry class\n");

  // Initialize the geometry/solver class for the two input cubes
  LrocPairModel lrocClass(params.leftFilePath, params.rightFilePath, params.worldTransform, params.includePosition);
  
  Vector<double> leftRow, leftCol, rightRow, rightCol;
  
  if (params.matchingPointsPath.empty())
  {
    printf("Searching for matching pixels in image overlap region\n");
    if (!findMatchingPixels(params, leftRow, leftCol, rightRow, rightCol))
      return false;
  }
  else
  {
    printf("Loading list of matched pixels from file\n");
    if (!loadMatchingPixels(params.matchingPointsPath, leftRow, leftCol, rightRow, rightCol))
      return false;
  }


  const size_t numMatchedPts = leftRow.size();

  // Load the inital points into the solver
  printf("Initializing solver state...\n");
  Vector<double> initialState, packedObservations;
  if (!lrocClass.getInitialStateEstimate(leftRow, leftCol, rightRow, rightCol, initialState, packedObservations, initialValues))
  {
    printf("Failed to get initial state!\n");
    return false;
  }


  //TODO: Move these
  std::string initialErrorPath    = "/home/smcmich1/initialAngleError.csv";
  std::string finalErrorPath      = "/home/smcmich1/finalAngleError.csv";
  std::string initialStatePath    = "/home/smcmich1/initialAngleState.csv";
  std::string finalStatePath      = "/home/smcmich1/finalAngleState.csv";
  std::string finalStateDiffPath  = "/home/smcmich1/finalAngleStateDiff.csv";
  std::string initialGdcCoordPath = "/home/smcmich1/initialGdcPoints.csv";
  std::string finalGdcCoordPath   = "/home/smcmich1/finalGdcPoints.csv";
  
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
   
  
  //printf("Writing initial state log to %s\n", initialStatePath.c_str());
  std::ofstream initialStateFile(initialStatePath.c_str()); 
  for (size_t i=0; i<initialState.size(); ++i)
    initialStateFile << initialState[i] << std::endl;
  initialStateFile.close();
  
  // Write initial points as GDC coordinates for google earth
  std::ofstream initialGdcCoordFile;
  if (params.initialOnly) // Write to output GDC file
    initialGdcCoordFile.open(params.gdcPointsOutPath.c_str());
  else // Use debug file path
    initialGdcCoordFile.open(initialGdcCoordPath.c_str());
  for (size_t i=0; i<numMatchedPts; ++i)
  {
    // Convert from GCC to GDC
    Vector3 gccPoint(initialState[startOfPts+ 3*i], initialState[startOfPts+ 3*i+1], initialState[startOfPts+ 3*i+2]);
    Vector3 gdcCoord = datum.cartesian_to_geodetic(gccPoint);
    initialGdcCoordFile.precision(12);
    // Write lat, lon, height so pc_align tool can read these files
    initialGdcCoordFile << gdcCoord[1] << ", " << gdcCoord[0] << ", " << gdcCoord[2] << std::endl;
  }
  initialGdcCoordFile.close();
  
  // Compute the initial error - euclidean point distance
  //printf("Writing initial error log to %s\n", initialErrorPath.c_str());
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
  printf("Mean point error before optimization = %lf\n", meanInitialError);
  
  Vector<double> initialComputedObservations = lrocClass(initialState);
  
  if (params.initialOnly)
  {
    printf("Stopping after initial state calculation\n");
    return true;
  }
  
   
  
  printf("Running solver...\n");
  Vector<double> finalParams;
  //@@@@@@@@@@@@@@@@@@@@ NEW METHOD @@@@@@@@@@@@@@@@@@@@@@@@@
#if false
  printf("Using Ceres solver!\n");
  //TODO: Turn on Glog
  
  // Create Ceres solver object
  ceres::Problem problem;
  
  // Set up camera parameters for solver
  double* cameraParams = &(initialState[0]);
  problem.AddParameterBlock(cameraParams, startOfPts);
  printf("startOfPts = %d\n", startOfPts);
  if (params.includePosition)
    printf("Six point cost function\n");
  
  for (size_t i=0; i<numMatchedPts; ++i) // For each input point
  {
    const int NUM_PARAMS_PER_POINT       = 3;
    const int NUM_PARAMS_PER_OBSERVATION = 2;
    
    // Set up this point's parameters for the solver
    double* pointParams = &(initialState[startOfPts + 3*i]);
    problem.AddParameterBlock(pointParams, NUM_PARAMS_PER_POINT);
    
    // Add the function and residual block for the left camera
    ceres::CostFunction* costFunctionLeft = 
            new ceres::NumericDiffCostFunction<LeftCostFunctor, ceres::CENTRAL, NUM_PARAMS_PER_OBSERVATION, 3>(
                new LeftCostFunctor(&lrocClass, leftCol[i], leftRow[i]));

    problem.AddResidualBlock(costFunctionLeft,
                             0 /* squared loss TODO experiment with this */,
                             pointParams);

    // Add the function and residual block for the right camera (slightly more complex
    ceres::CostFunction* costFunctionRight;
    if (params.includePosition) // Rotation and translation = 6 parameters
    {
      costFunctionRight = 
          new ceres::NumericDiffCostFunction<RightCostFunctor, ceres::CENTRAL, NUM_PARAMS_PER_OBSERVATION, 6, NUM_PARAMS_PER_POINT>(
                                new RightCostFunctor(&lrocClass, rightCol[i], rightRow[i]));
    }
    else // Translation only = 3 parameters
    {
      costFunctionRight = 
          new ceres::NumericDiffCostFunction<RightCostFunctor, ceres::CENTRAL, NUM_PARAMS_PER_OBSERVATION, 3, NUM_PARAMS_PER_POINT>(
                                new RightCostFunctor(&lrocClass, rightCol[i], rightRow[i]));
    }
    problem.AddResidualBlock(costFunctionRight,
                             0 /* squared loss TODO experiment with this */,
                             cameraParams, pointParams);
    
  } // End of loop through points
  
  // TODO: Select solver options
  ceres::Solver::Options solverOptions;
  solverOptions.linear_solver_type           = ceres::DENSE_SCHUR;
  solverOptions.minimizer_progress_to_stdout = true;
  
  // Execute the Ceres solver
  printf("Starting the Ceres solver...\n");
  ceres::Solver::Summary summary;
  ceres::Solve(solverOptions, &problem, &summary);
  
  std::ofstream ceresLog("~/data/ceresOutput.txt");
  std::cout << summary.FullReport() << "\n";
  ceresLog << summary.FullReport();
  ceresLog.close();

  printf("All done with the Ceres solver!\n");

  printf("Copying out Ceres solver results.\n");
  // Get the information back from the solver!
  finalParams.set_size(initialState.size());
  for (size_t i=0; i<initialState.size(); ++i)
    finalParams[i] = initialState[i];
  
  //return true;
  
#else
  //@@@@@@@@@@@@@@@@@@@@ OLD METHOD @@@@@@@@@@@@@@@@@@@@@@@@@
  printf("Using old solver!\n");
  // Now pass geometry class into the solver function
  int status;
  finalParams = vw::math::levenberg_marquardt(lrocClass, initialState, packedObservations, status); 
  
  std::cout << "Status = " << status << std::endl;
#endif
  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  
  //printf("Writing final error log to %s\n", finalErrorPath.c_str());
  //std::ofstream finalErrorFile(finalErrorPath.c_str()); 
  double meanFinalError = 0;
  std::vector<double> finalError = lrocClass.computeError(finalParams);
  Vector<double> finalPredictions = lrocClass(finalParams);
  Vector<double> rawError(finalPredictions.size());
  //std::ofstream predictionFile("/home/smcmich1/finalPredictions.csv");
  for (size_t i=0; i<finalPredictions.size(); ++i)
  {
  //  predictionFile << finalPredictions[i] << std::endl;
    rawError[i] = packedObservations[i] - finalPredictions[i];
    //finalErrorFile << packedObservations[i] - finalPredictions[i] << std::endl;
    if (i % 4 == 0)
      meanFinalError += finalError[i/4];
  }
  
  for (size_t i=0; i<finalError.size(); ++i)
  {   
  //  finalErrorFile << finalError[i] << std::endl;
    meanFinalError += finalError[i];
  }
  meanFinalError = meanFinalError / finalError.size(); 
  
  
  //predictionFile.close();
  //finalErrorFile.close();
  //meanFinalError = meanFinalError / currentError.size(); 
  
  //printf("Writing final state log to %s\n", finalStatePath.c_str());
  //std::ofstream finalStateFile(finalStatePath.c_str()); 
  //for (size_t i=0; i<finalParams.size(); ++i)
  //  finalStateFile << finalParams[i] << std::endl;
  //finalStateFile.close();

  // Write initial points as GDC coordinates for google earth
  std::ofstream finalGdcCoordFile;
  if (params.gdcPointsOutPath.empty())
    finalGdcCoordFile.open(finalGdcCoordPath.c_str()); // TODO: Remove debug default
  else
    finalGdcCoordFile.open(params.gdcPointsOutPath.c_str());
  for (size_t i=0; i<numMatchedPts; ++i)
  {
    // Convert from GCC to GDC
    Vector3 gccPoint(finalParams[startOfPts+ 3*i], finalParams[startOfPts+ 3*i+1], finalParams[startOfPts+ 3*i+2]);
    Vector3 gdcCoord = datum.cartesian_to_geodetic(gccPoint);
    finalGdcCoordFile.precision(12);
    // Write lat, lon, height so pc_align tool can read these files
    finalGdcCoordFile << gdcCoord[1] << ", " << gdcCoord[0] << ", " << gdcCoord[2] << std::endl;
  }
  finalGdcCoordFile.close();

  std::ofstream finalStateDiffFile(finalStateDiffPath.c_str()); 
  for (size_t i=0; i<finalParams.size(); ++i)
    finalStateDiffFile << finalParams[i] - initialState[i] << std::endl;
  finalStateDiffFile.close();
  
  printf("Mean point error after optimization = %lf\n", meanFinalError);
  printf("Mean error change = %lf\n", meanFinalError - meanInitialError);
  
  const double rad2deg = 180.0 / M_PI;
  const double deg2rad = M_PI / 180.0;
  printf("Adjustment rotation angles (degrees): %lf, %lf, %lf\n", finalParams[0]*rad2deg, finalParams[1]*rad2deg, finalParams[2]*rad2deg);
  //printf("Output rotation angles (degrees): %lf, %lf, %lf\n", finalParams[0]*rad2deg, finalParams[1]*rad2deg, finalParams[2]*rad2deg);
  if (params.includePosition)
    printf("Adjustment offset (meters): %lf, %lf, %lf\n", finalParams[3], finalParams[4], finalParams[5]);
  
  
  // These are instrument rotations from the frame file - rotation order is R = Rx(1.29)*Ry(-0.079)*Rz(180)
  // sc_from_instrument,  ( 1.129, -0.079, 180.0 ) [units in degrees]
  
  // vw::math::euler_to_rotation_matrix(a, b, c, 'order') --> R(c)*R(b)*R(a), order specifies which rotation axes
  
  //vw::math::Matrix<double,3,3> frameRotation = vw::math::euler_to_rotation_matrix(180.0*deg2rad, -0.079*deg2rad, 1.129*deg2rad, "zyx");
  vw::math::Matrix<double,3,3> frameRotation = vw::math::euler_to_rotation_matrix(1.129*deg2rad, -0.079*deg2rad, 180.0*deg2rad, "xyz");
   std:: cout << "Base rotation angles (degrees) = Rx(1.129) * Ry(-0.079) * Rz(180.0)"  << std::endl;

  // Existing offset matrix (from the fk file) = sc_from_instrument
  
  // This is the rotation matrix we solved for (matches code in IsisInterfaceLineScan.cc)
  // solvedRotation = instrument(good)_from_instrument
  vw::math::Matrix<double,3,3> solvedRotation = vw::math::euler_to_rotation_matrix(finalParams[0], finalParams[1], finalParams[2], "xyz");
  std::cout << "Solved rotation: " << solvedRotation << std::endl;
  // Multiply them to get the net rotation

  //vw::math::Matrix<double,3,3> outputRotation = frameRotation; //TEST: Make sure original rotations kept intact!
  //vw::math::Matrix<double,3,3> fixed_from_sc  = transpose(solvedRotation) * frameRotation;
  //vw::math::Matrix<double,3,3> outputRotation = transpose(fixed_from_sc); // == sc_from_fixed
  
  vw::math::Matrix<double,3,3> outputRotation;
  if (params.worldTransform)
    outputRotation = solvedRotation;
  else
    outputRotation = transpose(frameRotation * transpose(solvedRotation)); // == sc_from_instrument(good)
  std::cout << "Output rotation: " << outputRotation << std::endl;
  
  //vw::math::Matrix<double,3,3> outputRotationT = transpose(frameRotation * transpose(solvedRotation));
  //std::cout << "Output rotation TRANS: " << outputRotationT << std::endl;
  
  // Extract rotation angles --> Returns Rx(a)*Ry(b)*Rz(c)
  // -- Angle order needs to be  1, 2, 3
  //Vector3 outputAngles = vw::math::rotation_matrix_to_euler_zyx(outputRotation); //bad
  Vector3 outputAngles = vw::math::rotation_matrix_to_euler_xyz(outputRotation); // good?
  // Output angles (rZ, rY, rX)
  
  printf("Combined rotation angles (degrees): %lf, %lf, %lf\n", outputAngles[2]*rad2deg, outputAngles[1]*rad2deg, outputAngles[0]*rad2deg);
  
  if (!params.outputPath.empty()) // Dump rotation angles to a simple text file (in degrees)
  {
    printf("Writing output file %s\n", params.outputPath.c_str()); 
    std::ofstream outputFile(params.outputPath.c_str());
    //outputFile << outputAngles[0]*rad2deg << std::endl; // X rotation
    //outputFile << outputAngles[1]*rad2deg << std::endl; // Y rotation
    //outputFile << outputAngles[2]*rad2deg << std::endl; // Z rotation

    if (params.worldTransform) // Write the three (global) rotation angles as they are stored
    {
      for (int i=0; i<3; ++i)
        outputFile << finalParams[i] << std::endl;
    }
    else // Write out the combined relative rotation matrix
    {
      for (unsigned int r=0; r<outputRotation.rows(); ++r)
      {
        for (unsigned int c=0; c<outputRotation.cols(); ++c)
        {
          outputFile << outputRotation(r,c) << std::endl;
        }
      }
    }
    
    if (params.includePosition) // Also write out the translation
    {
      outputFile << finalParams[3] << std::endl;
      outputFile << finalParams[4] << std::endl;
      outputFile << finalParams[5];
    }
    outputFile.close();
    
    //TODO: Clean up this output section!
    if (params.worldTransform && params.includePosition)
    {
      // Write out the global rotation/translation matrix that was applied to the camera
      vw::math::Vector<double,3>   axisAngle   (finalParams[0], finalParams[1], finalParams[2]);
      vw::math::Quaternion<double> rotQuat     (axis_angle_to_quaternion(axisAngle));
      vw::math::Matrix<double,3,3> globalRotMat(rotQuat.rotation_matrix());
      
      std::string backupTransformFile(params.outputPath + ".matrix.csv");
      std::ofstream outputFileM(backupTransformFile.c_str());
      for (int q=0; q<3; ++q)
        outputFileM << globalRotMat[q][0] << " " << globalRotMat[q][1] << " " << globalRotMat[q][2] << " " << finalParams[3+q] << std::endl;
      outputFileM << "0 0 1 0" << std::endl;
      
      outputFileM.close();
    } // End world transform with position case
    
  } // End output path case
  

  return true;
}


//-------------------------------------------------------------------------------------------

int main(int argc, char* argv[]) 
{ 
  try 
  {
    printf("Using Ceres solver!\n");
    // Parse the input parameters
    Parameters params;
    if (!handle_arguments(argc, argv, params))
    {
      printf("Failed to parse input parameters!\n");
      return false;
    }
    
    optimizeRotations(params);
  } ASP_STANDARD_CATCHES;

  return 0;
}




