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

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <vw/InterestPoint/InterestData.h>


//#include <asp/Tools/stereo.h>
#include <stereo.h> // Using local copy

/// \file matchBinaryToCsv.cc Converts a binary .match file to an easy to read .csv file





int main( int argc, char *argv[] ) 
{
  std::string inputPath, outputPath="";
  
  po::options_description general_options("Options");
  general_options.add_options()
    ("help,h",        "Display this help message")  
    ("input,i",   po::value<std::string>(&inputPath),                "The input .match file")
    ("output,o",  po::value<std::string>(&outputPath)->default_value(""), "The output .csv file (default append .csv)");
    
  po::positional_options_description positional_desc;
  positional_desc.add("input",  1);
  
  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <input path>" << std::endl << std::endl;
  usage << general_options << std::endl;

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(general_options).positional(positional_desc).run(), vm );
    po::notify( vm );
  } catch (const po::error& e) {
    std::cout << "An error occured while parsing command line arguments.\n";
    std::cout << "\t" << e.what() << "\n\n";
    std::cout << usage.str();
    return 1;
  }


  try {
    
    // Read in the interest points from the binary file
    printf("Loading interest points\n");
    std::vector<vw::ip::InterestPoint> ip1, ip2;
    vw::ip::read_binary_match_file(inputPath, ip1, ip2);
    const size_t numPoints = ip1.size();
    if (ip2.size() != numPoints)
    {
      printf("Error - Got different numbers of interest points!\n");
      return 0;
    }
    
    // Open output file
    std::string outputPathReal = inputPath + ".csv";
    if (!outputPath.empty())
      outputPathReal = outputPath;
    std::ofstream outputFile(outputPathReal.c_str());
    if (outputFile.fail())
    {
      printf("Failed to open output file %s for writing\n", outputPathReal.c_str());
      return 0;
    }
    
    printf("Writing out interest points\n");
    for (size_t i=0; i<numPoints; ++i)
      outputFile << ip1[i].x << ", " << ip1[i].y << ", " << ip2[i].x << ", " << ip2[i].y << std::endl; 

    // Done writing the output file
    outputFile.close();

    printf("Wrote %d pairs of interest points\n", numPoints);
    
  }
  catch (const vw::Exception& e) 
  {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  return 0;
}










