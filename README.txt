
Build instructions:


Running the LRONAC Pipeline:



======= Overview of files =======

----- Configuration files -----

ceresInstallScript.sh = Notes to help in installing CERES on a Linux machine with no admin access
CMakeLists.txt = Main build script for lronacPipeline
FindStereoPipeline.cmake = Cmake function to help with the build.
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
lronacAngleSolver.cc        = Older version of the SBA tool using the CERES solver.  Fewer options but still in use.
lronacSolverModelDouble.h   = Camera model code for new SBA tool.
lronacSolverModel.h         = Camera model code for old SBA tool.
lronacSolverSupport.h       = Support code for both SBA tools.
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




