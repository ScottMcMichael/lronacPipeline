#!/bin/bash



# Making sure the old example still works
#/home/smcmich1/repo/StereoPipeline/src/asp/Tools/lronacAngleSolver --outputPath dummy.txt --gdcPointsOutPath /home/smcmich1/data/angleCorrectionTest/pairGdcCheck.csv /home/smcmich1/data/angleCorrectionTest/M123514622LE.corrected.cub /home/smcmich1/data/angleCorrectionTest/M123514622RE.corrected.cub

#/home/smcmich1/repo/StereoPipeline/src/asp/Tools/lronacAngleSolver --outputPath dummy.txt --gdcPointsOutPath /home/smcmich1/data/angleCorrectionTest/pairGdcCheck.csv /home/smcmich1/data/angleCorrectionTest/M123514622LE.lronaccal.lronacecho.cub /home/smcmich1/data/angleCorrectionTest/M123514622RE.lronaccal.lronacecho.cub

# --> Output can be similar without the position correction so a more detailed numerical comparison would be needed to tell the difference

# --> The solution is very sensitive to the set of points used to solve.  Unfortunately after 70 or so points the non-sparse solver becomes very slow.  Multiple trials or a better solver may be needed to improve results.

# --> The ground points from the angle solver also have a lot of variation, maybe very sensitive to pixel match quality?

# Output from M120168714 with translation only is really good!

#***********************************************************************************
# Testing the big stereo calibration program
~/repo/lronacPipeline/stereoCalibrationProcess.py --left /home/smcmich1/data/ARISTARCHU2/M120168714LE.IMG --right /home/smcmich1/data/ARISTARCHU2/M120168714RE.IMG --stereo-left /home/smcmich1/data/ARISTARCHU2/M120175500LE.IMG --stereo-right /home/smcmich1/data/ARISTARCHU2/M120175500RE.IMG --lola ~/data/ARISTARCHU2/RDR_310E312E_23N26NPointPerRow_csv_table.csv --keep --outputL ~/data/stereoCorrectionTest/M120168714LE.final.cub --outputR ~/data/stereoCorrectionTest/M120168714RE.final.cub --outputSL ~/data/stereoCorrectionTest/M120175500LE.final.cub --outputSR ~/data/stereoCorrectionTest/M120175500RE.final.cub --workDir ~/data/stereoCorrectionTest/workingDir

# Draw test output
~/repo/lronacPipeline/calibrationReport.py --input ~/data/stereoCorrectionTest/workingDir/pairGdcCheckPre.csv --output ~/data/stereoCorrectionTest/pairGdcCheckPreCeres700.kml --name pairCheckPreCeres700 --color blue

#~/repo/lronacPipeline/calibrationReport.py --input ~/data/stereoCorrectionTest/workingDir/pairGdcCheckInitial.csv --output ~/data/stereoCorrectionTest/pairGdcCheckInitialC.kml --name pairCheckInitialC --color blue

#~/repo/lronacPipeline/calibrationReport.py --input ~/data/stereoCorrectionTest/workingDir/pairGdcCheckFinal.csv --output ~/data/stereoCorrectionTest/pairGdcCheckFinalC.kml --name pairCheckFinalC --color blue

#~/repo/lronacPipeline/calibrationReport.py --input ~/data/stereoCorrectionTest/workingDir/gdcPointsCheckFinalRot.csv --output ~/data/stereoCorrectionTest/gdcPointsGdcCheckFinalRotC.kml --name gdcPointsCheckFinalRotC --color red

#~/repo/lronacPipeline/calibrationReport.py --input ~/data/stereoCorrectionTest/workingDir/gdcPointsCheckFinalRotStereo.csv --output ~/data/stereoCorrectionTest/gdcPointsGdcCheckFinalRotCStereo.kml --name gdcPointsCheckFinalRotCStereo --color red

# Draw input to pc_align
#~/repo/lronacPipeline/calibrationReport.py --input /home/smcmich1/data/stereoCorrectionTest/workingDir/gdcPointsLarge.csv --output ~/data/stereoCorrectionTest/inputGdcPointsL.kml --name inputGdcPointsL --color red --skip 100000

# Draw pc_align moved points
#~/repo/lronacPipeline/calibrationReport.py --input /home/smcmich1/data/stereoCorrectionTest/workingDir/pcAlignOutput/dem-trans_source.csv --output ~/data/stereoCorrectionTest/movedGdcPointsLT.kml --name movedGdcPointsLT --color blue --skip 10000

#~/repo/lronacPipeline/calibrationReport.py --input /home/smcmich1/data/stereoCorrectionTest/workingDir/pcAlignOutput/dem-trans_reference.csv --output ~/data/stereoCorrectionTest/movedGdcPointsL.kml --name movedGdcPointsL --color blue --skip 100000

#***********************************************************************************


# Generate comparison data
#mkdir /home/smcmich1/data/angleCorrectionTest/oldResults
#lronac2mosaic.py -t 4 -o /home/smcmich1/data/angleCorrectionTest/oldResults --keep ~/data/M112646261LE.IMG ~/data/M112646261RE.IMG
#exit 0


# Convert control point file into usable format
#/home/smcmich1/repo/lronacPipeline/extractQtieControlPoints.py --cnetPath ~/data/ARISTARCHU2/control_pointreg6.net --output controlPoints.csv



# Noproj the corrected data
#noproj from=/home/smcmich1/data/stereoCorrectionTest/M120168714LE.final.cub    to=/home/smcmich1/data/stereoCorrectionTest/M120168714LE.final.noproj.cub    match=/home/smcmich1/data/stereoCorrectionTest/M120168714LE.final.cub specs=/home/smcmich1/data/stereoCorrectionTest/noprojInstruments_fullRes.pvl

#noproj from=/home/smcmich1/data/stereoCorrectionTest/M120168714RE.final.cub    to=/home/smcmich1/data/stereoCorrectionTest/M120168714RE.final.noproj.cub    match=/home/smcmich1/data/stereoCorrectionTest/M120168714LE.final.cub specs=/home/smcmich1/data/stereoCorrectionTest/noprojInstruments_fullRes.pvl

#cp /home/smcmich1/data/stereoCorrectionTest/M120168714LE.final.noproj.cub /home/smcmich1/data/stereoCorrectionTest/M120168714LE.mosaic.cub

#handmos from=/home/smcmich1/data/stereoCorrectionTest/M120168714RE.final.noproj.cub mosaic=/home/smcmich1/data/stereoCorrectionTest/M120168714LE.mosaic.cub outsample=0 outline=0 matchbandbin=FALSE priority=ontop

#cubenorm from=/home/smcmich1/data/stereoCorrectionTest/M120168714LE.mosaic.cub to=/home/smcmich1/data/stereoCorrectionTest/M120168714LE.mosaic.norm.cub


# Do the same for the stereo pair
#noproj from=/home/smcmich1/data/stereoCorrectionTest/M120175500LE.final.cub    to=/home/smcmich1/data/stereoCorrectionTest/M120175500LE.final.noproj.cub    match=/home/smcmich1/data/stereoCorrectionTest/M120175500LE.final.cub specs=/home/smcmich1/data/stereoCorrectionTest/noprojInstruments_fullRes.pvl

#noproj from=/home/smcmich1/data/stereoCorrectionTest/M120175500RE.final.cub    to=/home/smcmich1/data/stereoCorrectionTest/M120175500RE.final.noproj.cub    match=/home/smcmich1/data/stereoCorrectionTest/M120175500LE.final.cub specs=/home/smcmich1/data/stereoCorrectionTest/noprojInstruments_fullRes.pvl

#cp /home/smcmich1/data/stereoCorrectionTest/M120175500LE.final.noproj.cub /home/smcmich1/data/stereoCorrectionTest/M120175500LE.mosaic.cub

#handmos from=/home/smcmich1/data/stereoCorrectionTest/M120175500RE.final.noproj.cub mosaic=/home/smcmich1/data/stereoCorrectionTest/M120175500LE.mosaic.cub outsample=0 outline=0 matchbandbin=FALSE priority=ontop

#cubenorm from=/home/smcmich1/data/stereoCorrectionTest/M120175500LE.mosaic.cub to=/home/smcmich1/data/stereoCorrectionTest/M120175500LE.mosaic.norm.cub



# Generate a DEM from the two mosaic images
#3parallel_stereo --corr-timeout 400 --alignment affineepipolar --subpixel-mode 1 --disable-fill-holes /home/smcmich1/data/stereoCorrectionTest/M120168714LE.mosaic.norm.cub /home/smcmich1/data/stereoCorrectionTest/M120175500LE.mosaic.norm.cub ~/data/stereoCorrectionTest/finalStereo/output --processes 8 --threads-multiprocess 4 --threads-singleprocess 32 --compute-error-vector
#--nodes-list PBS_NODEFILE --processes 4 --threads-multiprocess 16 --threads-singleprocess 32

# Find out the center latitude of the mosaic
#point2dem --errorimage -o /home/smcmich1/data/stereoCorrectionTest/finalStereo/output-DEM.tif /home/smcmich1/data/stereoCorrectionTest/finalStereo/output-PC.tif -r moon --tr 1 --t_srs "+proj=eqc +lat_ts=24.62 +lat_0=0 +a=1737400 +b=1737400 +units=m" --nodata -32767

# Create a hillshade image to check if the central errors are gone
#hillshade /home/smcmich1/data/stereoCorrectionTest/finalStereo/output-DEM.tif -o /home/smcmich1/data/stereoCorrectionTest/demHillshade.tif


# Call the (hacked) pipeline
#~/repo/lronacPipeline/lronacPipeline.py --input-folder ~/data/ARISTARCHU2




