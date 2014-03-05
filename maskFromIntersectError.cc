// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <vw/Image/ImageView.h>
#include <vw/Image/MaskViews.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Image/Filter.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/FileIO/DiskImageView.h>
#include <vw/Cartography/GeoReference.h>


using namespace vw;

/// \file maskFromIntersectError.cc Generates a good pixel mask of all pixels with
///                                 intersection error below a threshold.


/// Functor to evaluate if the sum squared value of a pixel is below a threshold
/// - In this program it is evaluating intersection error and generating a mask.
class SumSquaredThreshFunctor  : public ReturnFixedType<PixelGray<uint8> >
{
public:
  /// Construct with the threshold.
  SumSquaredThreshFunctor(float maxErrorSquared) : m_maxErrorSquared(maxErrorSquared) {}

  /// Return 255 if pixel value is below a threshold, 0 otherwise.
  inline PixelGray<uint8> operator() (PixelRGB<float> const &p) const
  {
    float errorSquared = p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
    if (errorSquared < m_maxErrorSquared)
      return PixelGray<uint8>(255);
    else    // Error above a threshold
      return PixelGray<uint8>(0);
  }

private:
  float m_maxErrorSquared; ///< Squared error pixels will be compared to.
};


int main( int argc, char *argv[] ) {

  std::string inputImagePath="", outputImagePath="";

  float errorThreshold;
  po::options_description general_options("Options");
  general_options.add_options()
    ("help,h",           "Display this help message")
    ("errorThreshold,e", po::value<float      >(&errorThreshold)->default_value(5), "Set value over which pixels are flagged as invalid");

  po::options_description positional("");
  positional.add_options()
    ("input-image",    po::value(&inputImagePath),  "Path to input  image file")
    ("output-image",   po::value(&outputImagePath), "Path to output image file");

  po::positional_options_description positional_desc;
  positional_desc.add("input-image",  1);
  positional_desc.add("output-image", 1);

  std::string usage("[options] <input-image> <output-image>\n");
  po::variables_map vm;
  try {
    po::options_description all_options;
    all_options.add(general_options).add(positional);

    po::store( po::command_line_parser( argc, argv ).options(all_options).positional(positional_desc).style( po::command_line_style::unix_style ).run(), vm );

    po::notify( vm );
  } catch (po::error const& e) {
    vw::vw_throw( vw::ArgumentErr() << "Error parsing input:\n"
                  << e.what() << "\n" << usage << general_options );
  }

  if ( !vm.count("input-image") || !vm.count("output-image"))
    vw_throw( vw::ArgumentErr() << "Requires <input-image> and <output-image> inputs in order to proceed.\n\n"
              << usage << general_options );

  float errorThresholdSquared = errorThreshold*errorThreshold;

  try {
  
    // Set up input from file
    vw::DiskImageView<PixelRGB<float> > inputImage(inputImagePath);


    // Set up thresholding function
    ImageViewRef<PixelGray<vw::uint8> > mask = vw::per_pixel_filter(inputImage,
                                                                    SumSquaredThreshFunctor(errorThresholdSquared));


    // Set up output file
    boost::scoped_ptr<DiskImageResource> r(DiskImageResource::create(outputImagePath, mask.format()));

    // Copy georeference information from the input image
    vw::cartography::GeoReference georef;
    vw::cartography::read_georeference(georef, inputImagePath);
    write_georeference( *r, georef );

    // Do everything!
    vw::block_write_image( *r, mask,
                          vw::TerminalProgressCallback( "maskFromIntersectError", "Generating mask:") );

  }
  catch (const Exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  return 0;
}










