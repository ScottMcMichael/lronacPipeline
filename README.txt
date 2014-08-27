

======= Build instructions ===========================================

- Make sure your StereoPipeline installation is up to date.
  - Last tested with the ASP version listed blow:
"
NASA Ames Stereo Pipeline 2.4.0_post
  Build ID: ffbcda1
Built against:
  NASA Vision Workbench 2.2.0_post
    Build ID: 32b9d80
  USGS ISIS 3.4.6.5819 2014-04-29 # Version date v002 # 3rd party libraries version stable # release stage (alpha, beta, stable)
  Boost C++ Libraries 105500
  GDAL 1.10.1 | 20130826
  Proj.4 480
"

- Check out and build the NGT Tools repository:
git clone https://github.com/NeoGeographyToolkit/Tools.git
[ Build according to the instructions in the Tools/README.txt file. Make sure not to build in debug mode! ]

- Add NGT Tools to your PATH and PYTHONPATH variables
- Add NGT Tools build directory to PATH

- Check out the lronacPipeline repository:
git clone https://github.com/ScottMcMichael/lronacPipeline.git
cd lronacPipeline

- Edit CMakeLists.txt to point to the correct BaseSystem, visionworkbench, StereoPipeline, and NGT Tools installations.

- Build lronacPipeline:
mkdir build
cd build
cmake ..
make

- Add lronacPipeline and lronacPipeline/build to your path.



======= Running the LRONAC Pipeline =======================================

- To generate a fully corrected DEM, call lronac2dem.py similar to this:

lronac2dem.py --left /byss/moon/lronacPipeline_V2/BHABHAPLAIN/M112646261LE.IMG --right /byss/moon/lronacPipeline_V2/BHABHAPLAIN/M112646261RE.IMG --stereo-left /byss/moon/lronacPipeline_V2/BHABHAPLAIN/M112653051LE.IMG --stereo-right /byss/moon/lronacPipeline_V2/BHABHAPLAIN/M112653051RE.IMG --lola /byss/moon/lronacPipeline_V2/BHABHAPLAIN/lolaRdrPts.csv --workDir /byss/moon/lronacPipeline_V2/BHABHAPLAIN/workdir --outputPrefix /byss/moon/lronacPipeline_V2/BHABHAPLAIN/results/output

Note that you need the two pairs of LRONAC IMG files plus the LOLA data covering the expected region.  You can obtain the files manually or you can use lronacDataGrabber to find them on the ASU LRONAC website.



======= Downloading data ==================================================

- Use the lronacDataGrabber.py tool to generate a list of data locations.

lronacDataGrabber.py --fetch --input-file lronacDataSourceList.txt

- Then use the same script to retrieve data files from the list.

lronacDataGrabber.py --input-file lronacDataSourceList.txt --output-folder . --name BHABHAPLAIN

The minimum bounds for each data set are listed in the data grabber source location file and in the command line output when retrieving data.  A 0.25 degree buffer in each direction is reccomended.  


======= Output file description ====================================================

-ASU_LOLA_diff_stats.txt  = Contains statistics comparing elevation differences between the LOLA DEM and the ASU DEM.
-ASU_LOLA_diff_points.txt = List of each point compared to generate the statistics in the previous file (lat, lon, elevation, error).
-LOLA_diff_stats.txt      = Contains statistics comparing elevation differences between the LOLA DEM and the newly generated DEM.
-LOLA_diff_points.txt     = List of each point compared to generate the statistics in the previous file (lat, lon, elevation, error).
-Hillshade.tif        = Convenient 8-bit visualization of the output DEM.
-Colormap.tif         = Elevation colorized version of Hillshade.tif
-DEM.tif              = The output 32 bit floating point DEM.
-IntersectionErr.tif  = Image containing the stereo triangulation error of each pixel in the output DEM.
-MapProjLeft.tif      = Map projection of the main input files.
-MapProjLeft.tif      = Map projection of the stereo input files.

output-Log.txt = Log of intermediate processing outputs.  There are a number of things to look for in this file:
- For the two left cubes: "INFO:root:- Stereo completed with good pixel percentage: 0.843245660555".
    A low good pixel percentage means that the stereo algorithm may be having trouble with the input images.
- For the the other cube pairs: "INFO:root:- Found 653 point pairs"
    A low number means the interest point finding algorithm may be having trouble with the input images.
    For the cropped images 50 to 100 points is ok.  For full size images > 300 is best.
- The median point error before and after SBA computation (look for ">>>>").
    The final output error should be < 1.0 pixels.  A higher value means the SBA solver had trouble finding a solution.
- The mean pair errors (look for "=====>").
    These errors are computed using "campt" on pixel pairs of the fully corrected cube files at the end of the calibration process.
    The main and stereo pairs should have a small error (< 1.0 meters) but the other two errors could be higher (< 60 meters).
- The jitreg offsets ("INFO:root:- jitreg offsets = [0.086819999999999994, 0.97209999999999996]")
    These values should be < 2.0 pixels.
- There are timing results scattered throughout the file.  Some time values include several of the previous times.

======= Processing description =====================================================

--- Overview

1 - Apply known position offset (about 1 meter per camera)
2 - Run stereo program to find pixel matches between all camera pairs
3 - Put all pixel pairs through SBA solver to compute transforms between cameras
4 - Use pc_align to find offset of transformed cameras from the LOLA data.
5 - Apply all internal (SBA) and external (pc_align) transforms to cameras
6 - Generate mosaics from corrected cameras.
7 - Generate final stereo DEM from mosaics.


--- Detailed steps

Inputs: Lola point cloud, four images (LE1, RE1, LE2, RE2)
1  - Preprocess data with lronac2isis, lronacecho, spiceinit
2  - Apply known position offset (about 1 meter per camera)
3  - Run stereo program to find pixel matches between all camera pairs
4  - Put all pixel pairs through SBA solver to compute transforms between cameras
   - Solve for global rotation (RE1, LE2, RE2) and translation (2 from 1)
5  - Apply transforms to the cameras (RE1, LE2, and RE2)
6  - Generate a point cloud from the camera pair with the greatest overlap
7  - Use pc_align to find offset of the new point cloud to the LOLA point cloud
8  - Apply pc_align computed transform to all four cameras
9  - Noproj the two pairs of cameras
10 - Make a mosaic with each of the two pairs using lronacjitreg and handmos
11 - Generate final stereo DEM from mosaics


Note that the processing method changed after version 1.4 of the software.  The local rotation 
computations and generation of modified kernel files were removed to cope with the significant 
changes to the LRO-NAC kernel structure that came with the release of the newer temperature
dependent kernels.  Most of the 'V2' files were introduced to handle these changes.


======= Overview of files =====================================================

----- Configuration files -----

CMakeLists.txt = Main build script for lronacPipeline


----- Executables -----

--> Neither of these tools are used by the current version.
mkspk  = NAIF provided tool to create new position kernel files.
msopck = NAIF provided tool to create new orientation kernel files.


----- Python Scripts -----

rotationCorrector.py   = Calls mkspk and msopck to generate modified position and orientation kernel files. 
positionCorrector.py   = Calls mkspk to generate a modified position kernel file (applies rigid camera offset). 
rotationCorrectorV2.py = Modifies the position and orientation of a cube file.
positionCorrectorV2.py = Modifies the position of a cube file (applies rigid camera offset).


lronacCameraRotationCorrector.py = Applies a rotation to a frame kernel file.

lronacPipeline.py    = Old test script based on lronac2mosaic.py
lronac2dem.py        = Generates a calibrated DEM and packages it up for archiving.
IsisTools.py         = Python functions used by other scripts in this folder.
makeDemAndCompare.py = Generates a DEM and map projections from two cubes and calculates accuracy statistics.

stereoDoubleCalibrationProcess.py = Given two pairs of .IMG files, generates fully calibrated version of each of them.

lronacDataGrabber.py = Tool to help download batches of LRONAC data matching the ASU DTM data.
generateReport.py    = Tool for plotting the output of lronacPipeline.py.
calibrationReport.py = Generates a KML plot of a set of GDC points. 
archiveManager.py    = Specialized tool for handling data on the Nasa Pleiades computer.
clearOldKernels.py   = Script for deleting deprecated LRO-NAC spice kernels from an ISIS installation to save space.


----- C++ code -----

lronacAngleDoubleSolver.cc  = New SBA tool using the CERES solver.
lronacSolverModelDouble.h   = Camera model code for new SBA tool.
lronacSolverSupport.h       = Support code for the SBA tool.
IsisInterfaceLineScanRot.h  = Replacement of IsisInterfaceLineScan with additional functionality.
IsisInterfaceLineScanRot.cc = See IsisInterfaceLineScanRot.h.
SpiceEditor.cc              = Tool to generate intermediate position and orientation kernel correction files given a transform.
SpiceEditorV2.cc            = Tool to apply a correction to a cube's position/rotation data given a transform.
stereo.h                    = Copy of the same file from StereoPipeline/src/asp/Tools with some dependencies removed.
stereo.cc                   = See stereo.h

----- Unused/Debug files -----

makePcAlignPlots.py         = Tool to help visualize pc_align output.  Not currently used.
lronacMathTester.py         = Debug script.

