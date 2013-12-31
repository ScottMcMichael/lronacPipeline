

======= Build instructions ===========================================

- Make sure your StereoPipeline installation is up to date.

- Install the CERES solver library.  Look in ceresInstallScript.sh for hints if you need to install manually.

- Check out the lronacPipeline repository:
git clone https://github.com/ScottMcMichael/lronacPipeline.git
cd lronacPipeline

- Edit CMakeLists.txt to point to the correct BaseSystem, visionworkbench, StereoPipeline, CERES, and SuiteSparse installations.

- Build lronacPipeline:
mkdir build
cd build
cmake ..
make

- Add lronacPipeline and lronacPipeline/build to your path.



======= Running the LRONAC Pipeline =======================================

- To generate a fully corrected DEM, call lronac2dem.py similar to this:

lronac2dem.py --left /byss/moon/lronacPipeline_V2/BHABHAPLAIN/M112646261LE.IMG --right /byss/moon/lronacPipeline_V2/BHABHAPLAIN/M112646261RE.IMG --stereo-left /byss/moon/lronacPipeline_V2/BHABHAPLAIN/M112653051LE.IMG --stereo-right /byss/moon/lronacPipeline_V2/BHABHAPLAIN/M112653051RE.IMG --lola /byss/moon/lronacPipeline_V2/BHABHAPLAIN/lolaRdrPts.csv --workDir /byss/moon/lronacPipeline_V2/BHABHAPLAIN/workdir --prefix /byss/moon/lronacPipeline_V2/BHABHAPLAIN/results/output

Note that you need the two pairs of LRONAC IMG files plus the LOLA data covering the expected region.  You can obtain the files manually or you can use lronacDataGrabber to find them on the ASU LRONAC website.



======= Downloading data ==================================================

- Use the lronacDataGrabber.py tool to generate a list of data locations.

lronacDataGrabber.py --fetch --input-file lronacDataSourceList.txt

- Then use the same script to retrieve data files from the list.

lronacDataGrabber.py --input-file lronacDataSourceList.txt --output-folder . --name BHABHAPLAIN

Currently the server for automated LOLA data downloads is not working so they needed to be downloaded by hand at the following address:  http://ode.rsl.wustl.edu/moon/lrololadataPointSearch.aspx
The minimum bounds for each data set are listed in the data grabber source location file and in the command line output when retrieving data.  A 0.25 to 0.5 degree buffer in each direction is reccomended.  



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
   - Solve for local rotation (RE to LE) and global rotation (2 to 1)
5  - Apply global transform to the stereo camera pair (LE2 and RE2).
6  - Generate a point cloud from LE1 and LE2
7  - Use pc_align to find offset of the new point cloud to the LOLA point cloud
8  - Apply pc_align computed transform to all four cameras.
9  - Apply local transform to RE1 and RE2
10 - Noproj the two pairs of cameras
11 - Make a mosaic with each of the two pairs using lronacjitreg and handmos
12 - Generate final stereo DEM from mosaics.



======= Overview of files =====================================================

----- Configuration files -----

ceresInstallScript.sh     = Notes to help in installing CERES on a Linux machine with no admin access
CMakeLists.txt            = Main build script for lronacPipeline
FindStereoPipeline.cmake  = Cmake function to help with the build.
FindVisionWorkbench.cmake = Cmake function to help with the build.


----- Executables -----

mkspk  = NAIF provided tool to create new position kernel files.
msopck = NAIF provided tool to create new orientation kernel files.


----- Python Scripts -----

rotationCorrector.py = Calls mkspk and msopck to generate modified position and orientation kernel files. 
positionCorrector.py = Calls mkspk to generate a modified position kernel file. 

lronacCameraRotationCorrector.py = Applies a rotation to a frame kernel file.

lronacPipeline.py = Original test script based on lronac2mosaic.py
lronac2dem.py     = New tool to generate a more accurate DEM from two pairs of IMG files.
IsisTools.py      = Python functions used by other scripts in this folder.

stereoDoubleCalibrationProcess.py = Given two pairs of .IMG files, generates fully calibrated version of each of them.

lronacDataGrabber.py = Tool to help download batches of LRONAC data matching the ASU DTM data.
generateReport.py    = Tool for plotting the output of lronacPipeline.py.
calibrationReport.py = Generates a KML plot of a set of GDC points. 



----- C++ code -----

lronacAngleDoubleSolver.cc  = New SBA tool using the CERES solver.
lronacSolverModelDouble.h   = Camera model code for new SBA tool.
lronacSolverSupport.h       = Support code for the SBA tool.
pixelPairsFromStereo.cc     = Tool to extract a grid of correspondence points from a stereo output file.
IsisInterfaceLineScanRot.h  = Replacement of IsisInterfaceLineScan with additional functionality.
IsisInterfaceLineScanRot.cc = See IsisInterfaceLineScanRot.h.
SpiceEditor.cc              = Tool to generate intermediate position and orientation kernel correction files given a transform.
stereo.h                    = Copy of the same file from StereoPipeline/src/asp/Tools with some dependencies removed.
stereo.cc                   = See stereo.h
lola_compare.cc             = Tool for comparing a GeoTiff DEM to a LOLA RDR file.  Also found in NeoGeographyToolkit/Tools.
imagestats.cc               = Tool similar to gdalinfo -stats but provides better histogram output.  Also found in NeoGeographyToolkit/Tools.

----- Unused/Debug files -----

makePcAlignPlots.py         = Tool to help visualize pc_align output.  Not currently used.
lronacMathTester.py         = Debug script.
extractQtieControlPoints.py = Converts a Qtie control points file to a simple csv format.


