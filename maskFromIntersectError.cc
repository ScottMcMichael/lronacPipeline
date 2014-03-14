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
  /// Construct with the thresholds.
  SumSquaredThreshFunctor(const std::vector<float> &thresholdValues, bool scale=false) : m_scaleVal(1.0)
  {
    // Record squared copies of each threshold value
    m_ValsSquared.resize(thresholdValues.size());
    for (size_t i=0; i<m_ValsSquared.size(); ++i)
      m_ValsSquared[i] = thresholdValues[i]*thresholdValues[i];


    // Compute scale value
    // - Number of output levels = 1 + number of thresholds
    // - Outputs will be 0, 1*m_scaleVal, 2*m_scaleVal, ...
    if (scale)
      m_scaleVal = 255.0f / float(thresholdValues.size());
  }

  /// Assign each pixel a value based on the threshold values.
  inline PixelGray<uint8> operator() (PixelRGB<float> const &p) const
  {
    float valueSquared = p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
    uint8 bin = 0;
    for (size_t i=0; i<m_ValsSquared.size(); ++i)
    {
      if (valueSquared >= m_ValsSquared[i]) // Value exceeds this threshold
        bin = i+1;  // Set it to threshold index
      else  // Get out of the loop
        break;
    }
    // Now apply the scale value to the output
    return PixelGray<uint8>(bin*m_scaleVal);
  }

private:
  std::vector<float> m_ValsSquared;
  float m_scaleVal;
};


int main( int argc, char *argv[] ) {

  std::string inputImagePath="", outputImagePath="";

  std::vector<float> thresholdValues;
  bool scaleOutput;
  po::options_description general_options("Options");
  general_options.add_options()
    ("help,h",         "Display this help message")
    ("scaleOutput,s",  po::bool_switch(&scaleOutput)->default_value(false),  "Scale output up to 255")
    ("thresholds",     po::value<std::vector<float> >(&thresholdValues)->multitoken(),
                         "Set one or more pixel threshold values (must be >0 and increasing value)");

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

  if ( !vm.count("input-image") || !vm.count("output-image") || !vm.count("thresholds"))
    vw_throw( vw::ArgumentErr() << "Requires <input-image>, <output-image>, and threshold inputs in order to proceed.\n\n"
              << usage << general_options );

  // Verify that the threshold values are sorted and > 0
  float lastVal = 0;
  for (size_t i=0; i<thresholdValues.size(); ++i)
  {
    if (thresholdValues[i] <= lastVal)
      vw_throw( vw::ArgumentErr() << "Threshold values must be greater than zero and in increasing order.\n\n");
    lastVal = thresholdValues[i];
  }

  try {
  
    // Set up input from file
    vw::DiskImageView<PixelRGB<float> > inputImage(inputImagePath);


    // Set up thresholding function
    ImageViewRef<PixelGray<vw::uint8> > mask = vw::per_pixel_filter(inputImage,
                                                                    SumSquaredThreshFunctor(thresholdValues,
                                                                                            scaleOutput));

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










