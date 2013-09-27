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


/// \file lrojitreg.cc
///

#include <asp/Tools/stereo.h>
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
#include <asp/IsisIo/IsisCameraModel.h> //::point_to_pixel>
#include <vw/Math/LevenbergMarquardt.h>


using namespace vw;
using namespace vw::stereo;
using namespace asp;
using std::endl;
using std::setprecision;
using std::setw;


struct Parameters : asp::BaseOptions 
{
  // Input paths
  std::string leftFilePath;
  std::string rightFilePath;
  //std::string rowLogFilePath;

  // Settings
  //float log;
  //int   h_corr_min, h_corr_max;
  //int   v_corr_min, v_corr_max;
  //Vector2i kernel;
  //int   lrthresh;
  //int   correlator_type;
  int   cropWidth;  
};



bool handle_arguments(int argc, char* argv[],
                     Parameters &opt) 
{ 
  po::options_description general_options("Options");
  general_options.add_options()
    ("crop-width",       po::value(&opt.cropWidth )->default_value(200), "Crop images to this width before disparity search");
  
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

//----------------------------------------------------------------------------------------------------------------
/// Class for solving for the rotation between two LRONAC cameras
class LrocPairModel : public vw::math::LeastSquaresModelBase<LrocPairModel>
{
public:  // Definitions
  
  /// * Defines a result_type that is the type returned by
  ///   evaluating the functor.  Typically Vector<float> or
  ///   Vector<double>
  typedef result_type Vector<double>;

  /// * Defines a domain_type that is the type of the search
  ///   space.  Often a Vector<foo>, but can reflect other
  ///   topologies if needed.
  typedef domain_type Vector<double>;
  
  
  /// * Defines a jacobian_type corresponding to the space of
  ///   jacobian matrices.  Typically Matrix<foo>.
  typedef jacobian_type Matrix<double>;
  
private: // Variables
  
  // Camera models
  vw::camera::IsisCameraModel     _leftCameraModel;
  vw::camera::IsisCameraModel     _rightCameraModel;
  
  boost::shared_pointer<vw::camera::IsisCameraModel> _rightCameraPointer;
  vw::camera::AdjustedCameraModel _rightCameraRotatedModel;
  
  vw::stereo::StereoModel         _stereoModel;
  
  // Observation records
  Vector<double> _leftRows;
  Vector<double> _rightRows;
  Vector<double> _leftCols;
  Vector<double> _rightCols;
  
public: // Functions  
  
  /// Constructor performs initialization
  LrocPairModel(const std::string &leftCubePath, const std::string &rightCubePath)
    : _leftCameraModel(leftCubePath), _rightCameraModel(rightCubePath), 
      _rightCameraPointer(&rightCubePath), // TODO: Existing null deleter function?
      _rightCameraRotatedModel(_rightCameraPointer),
      _stereoModel(&_leftCameraModel, &_rightCameraRotatedModel) // Only the right camera is rotated
  {
    // Both camera models are loaded from file on initialization
  }
  
  
  /// Given the pixel pair observations, compute the initial state estimate.
  /// - This also loads the observation vectors and returns a packed version of them.
  bool getInitialStateEstimate(const Vector<double> &leftRows,  const Vector<double> &leftCols,
                               const Vector<double> &rightRows, const Vector<double> &rightCols,
                               Vector<double> &stateEstimate, Vector<double> &packedObsVector)
  { 
    //TODO: Verify vector sizes match
    
    // Record observation vectors to class variables
    _leftRows  = leftRows;
    _rightRows = rightRows;
    _leftCols  = leftCols;
    _rightCols = rightCols;    // Maybe don't need these

    /* Input parameter set:
      rotationOffsetX (additional rotations applied to RE camera)
      rotationOffsetY
      rotationOffsetZ
      x               [Repeated for every correspondence point]
      y
      z
    */
    
    // Determine the number of state elements
    const size_t numPoints        = leftRows.size();
    const size_t numStateElements = numPoints*3 + 3;
    stateEstimate.resize(numStateElements);
    packedObsVector.resize(numPoints*4);
    
    // Rotation values start at zero
    stateEstimate[0] = 0;
    stateEstimate[1] = 0;
    stateEstimate[2] = 0;

    // For each input point pair
    Vector3 pointLoc;
    for (size_t i=0; i<numPoints; ++i)
    {
      // Set up point pair
      Vector2<double> leftPixel (leftCols[i],  leftRows[i]);
      Vector2<double> rightPixel(rightCols[i], rightRows[i]);
      
      // Compute the intersection location
      double triangulationError;
      pointLoc = _stereoModel(leftPixel, rightPixel, triangulationError);
      
      // Record the x/y/z value for this point
      stateEstimate[3 + i*3 + 0] = pointLoc[0];
      stateEstimate[3 + i*3 + 1] = pointLoc[1];
      stateEstimate[3 + i*3 + 2] = pointLoc[2];
      
      // Build the packed observation vector
      packedObsVector[i*4 + 0] = leftCols [i];
      packedObsVector[i*4 + 1] = leftRows [i];
      packedObsVector[i*4 + 2] = rightCols[i];
      packedObsVector[i*4 + 3] = rightRows[i];
    } // End loop through points
    
    return true;
  } // end getInitialStateEstimate()
  
  /// * Defines a method: result_type operator()( domain_type const& x ) const;
  ///   that evaluates the model function at the given point.
  result_type operator()( domain_type const& x ) const
  {
    // This function returns an error vector for a given set of parameters

    //TODO: Verify all rotation orders etc!
    // Apply the rotations from the state vector to the right LROC camera model
    _rightCameraRotatedModel.set_axis_angle_rotation(Vector3<double>(x[0], x[1], x[2]));
    
    // Set up output vector
    size_t numPoints = _leftRows.size();
    Vector<double> obsVec(numPoints*4);
    
    // Compute expected obseration value at each pixel
    for (size_t i=0; i<numPoints; ++i)
    {
      // Create a point object
      Vector3<double> thisPoint(x[3 + observationNum*3 + 0],
                                x[3 + observationNum*3 + 1]
                                x[3 + observationNum*3 + 2]);
                                
      // Project the point into both cameras
      Vector2<double> leftProjection  = _leftCameraModel.point_to_pixel(thisPoint);
      Vector2<double> rightProjection = _rightCameraRotatedModel.point_to_pixel(thisPoint);
      
      // Load the projected pixels into the output obseration vector
      obsVec[4*i + 0] = leftProjection [0]; // x
      obsVec[4*i + 1] = leftProjection [1]; // y
      obsVec[4*i + 2] = rightProjection[0];
      obsVec[4*i + 3] = rightProjection[1];
    }
    
    return obsVec;
  }
  
  
  /// * The domain_type must implement a method: domain_type domain_type::operator+( gradient_type const& g ) const;
  ///   that adds a tangent vector to a position.  You get this for
  ///   free if both domain_type and gradient_type are Vector<foo>.
  ///   This is where you do most of the hard work if domain_type
  ///   represents some non-trivial topological space.
  
  
  /// * The result_type must implement a method: double result_type::norm_2( result_type const& g ) const;
  ///   that is used in some optimizers to compute the error.  You get this for
  ///   free if result_type is a Vector<foo>.
  
  
  /// * The jacobian_type must implement several matrix-like methods such as
  ///   scalar multiplication on the left.  You get this for free in the usual case when
  ///   jacobian_type is just Matrix<foo>.
  ///
  
}; // End class LrocPairModel
//----------------------------------------------------------------------------------------------------------------

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
  
  
  // Load both images  
  printf("Loading images left=%s and right=%s...\n",
         params.leftFilePath.c_str(),
         params.rightFilePath.c_str());
  DiskImageView<PixelGray<float> > left_disk_image (params.leftFilePath );
  DiskImageView<PixelGray<float> > right_disk_image(params.rightFilePath);
  
  
  const int imageWidth      = std::min(left_disk_image.cols(), right_disk_image.cols());
  const int imageHeight     = std::min(left_disk_image.rows(), right_disk_image.rows());
  const int imageTopRow     = 0;
  const int imageMidPointX  = imageWidth / 2;
  const int cropStartX      = imageMidPointX - (params.cropWidth/2);

  // Restrict processing to the border of the images
  // - Since both images were nproj'd the overlap areas should be in about the same spots.
  const BBox2i crop_roi( cropStartX, imageTopRow,
			 params.cropWidth, imageHeight );
  std::cout << "Expected overlap ROI = " << crop_roi << std::endl;

  const int SEARCH_RANGE_EXPANSION = 5;
  int  ipFindXOffset = 0;
  int  ipFindYOffset = 0;
  bool ransacSuccess = false;

  // Now use interest point finding/matching functions to estimate the search offset between the images
  printf("Gathering interest points...\n");

  // Gather interest points
  asp::IntegralAutoGainDetector detector( 500 );
  ip::InterestPointList ip1 = ip::detect_interest_points( vw::create_mask_less_or_equal(crop(left_disk_image, crop_roi), 0), detector );
  ip::InterestPointList ip2 = ip::detect_interest_points( vw::create_mask_less_or_equal(crop(right_disk_image,crop_roi), 0), detector );
  printf("Found %lu, %lu interest points.\n", ip1.size(), ip2.size());
      
  ip::SGradDescriptorGenerator descriptor;
  describe_interest_points( vw::create_mask_less_or_equal(crop(left_disk_image, crop_roi), 0), descriptor, ip1 );
  describe_interest_points( vw::create_mask_less_or_equal(crop(right_disk_image,crop_roi), 0), descriptor, ip2 );

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
  
  printf("Found %lu, %lu matched interest points.\n", matched_ip1.size(), matched_ip2.size());

  
  
/* TODO: Retain only RANSAC matching points?
  // Filter interest point matches
  math::RandomSampleConsensus<math::SimilarityFittingFunctor, math::InterestPointErrorMetric> ransac( math::SimilarityFittingFunctor(),
                          math::InterestPointErrorMetric(),
                          100, 5, 100, true );
  std::vector<Vector3> ransac_ip1 = ip::iplist_to_vectorlist(matched_ip1);
  std::vector<Vector3> ransac_ip2 = ip::iplist_to_vectorlist(matched_ip2);

  // Finding offset using RANSAC...
  try
  {
    Matrix<double> H(ransac(ransac_ip1, ransac_ip2));
    std::cout << "ipfind based similarity: " << H << std::endl;

    // Use the estimated transform between the images to determine a search offset range
    ipFindXOffset = static_cast<int>(H[0][2]);
    ipFindYOffset = static_cast<int>(H[1][2]);
    ransacSuccess = true;
  }
  catch(...) // Handle a RANSAC failure
  {
    printf("RANSAC solution failed!\n");
    return false;
  }
*/

  // Now that we have correspondence points, feed them into an angular solver.

  // TODO: Convert the matching points into the correct format!
  Vector<double> leftRow, leftCol, rightRow, rightCol;

  printf("Constructing geometry class\n");

  // Initialize the geometry/solver class for the two input cubes
  LrocPairModel lrocClass(leftCube, rightCube);
  Vector<double> initialState, packedObservations;
  lrocClass.getInitialStateEstimate(leftRow, leftCol, rightRow, rightCol, initialState, packedObservations);

  printf("Running solver...\n");

  // Now pass geometry class into the solver function
  int status;
  Vector<double> finalParams = vw::math::levenberg_marquardt(lrocClass, initialState, packedObservations, status);
//                                                   double abs_tolerance = VW_MATH_LM_ABS_TOL,
//                                                   double rel_tolerance = VW_MATH_LM_REL_TOL,
//                                                   double max_iterations = VW_MATH_LM_MAX_ITER)

  printf("Output rotation angles: %lf, %lf, %lf\n", finalParams[0], finalParams[1], finalParams[2]);



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
    
    optimizeRotations(params);
  } ASP_STANDARD_CATCHES;

  return 0;
}






