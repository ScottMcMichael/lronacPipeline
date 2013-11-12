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
#~/repo/lronacPipeline/stereoCalibrationProcess.py --left /home/smcmich1/data/ARISTARCHU2/M120168714LE.IMG --right /home/smcmich1/data/ARISTARCHU2/M120168714RE.IMG --stereo-left /home/smcmich1/data/ARISTARCHU2/M120175500LE.IMG --stereo-right /home/smcmich1/data/ARISTARCHU2/M120175500RE.IMG --lola ~/data/ARISTARCHU2/RDR_310E312E_23N26NPointPerRow_csv_table.csv --keep --outputL ~/data/stereoCorrectionTest/M120168714LE.final.cub --outputR ~/data/stereoCorrectionTest/M120168714RE.final.cub --workDir ~/data/stereoCorrectionTest/workingDir

# Draw test output
#~/repo/lronacPipeline/calibrationReport.py --input ~/data/stereoCorrectionTest/workingDir/pairGdcCheckInitial.csv --output ~/data/stereoCorrectionTest/pairGdcCheckInitialC.kml --name pairCheckInitialC --color blue

#~/repo/lronacPipeline/calibrationReport.py --input ~/data/stereoCorrectionTest/workingDir/pairGdcCheckFinal.csv --output ~/data/stereoCorrectionTest/pairGdcCheckFinalC.kml --name pairCheckFinalC --color blue

#~/repo/lronacPipeline/calibrationReport.py --input ~/data/stereoCorrectionTest/workingDir/gdcPointsCheckFinalRot.csv --output ~/data/stereoCorrectionTest/gdcPointsGdcCheckFinalRotC.kml --name gdcPointsCheckFinalRotC --color red


# Draw input to pc_align
#~/repo/lronacPipeline/calibrationReport.py --input /home/smcmich1/data/stereoCorrectionTest/workingDir/gdcPointsLarge.csv --output ~/data/stereoCorrectionTest/inputGdcPoints.kml --name inputGdcPoints --color red --skip 100

# Draw pc_align moved points
#~/repo/lronacPipeline/calibrationReport.py --input /home/smcmich1/data/stereoCorrectionTest/workingDir/pcAlignOutput/dem-trans_source.csv --output ~/data/stereoCorrectionTest/movedGdcPoints.kml --name movedGdcPoints --color red --skip 100

#***********************************************************************************



# Generate comparison data
#mkdir /home/smcmich1/data/angleCorrectionTest/oldResults
#lronac2mosaic.py -t 4 -o /home/smcmich1/data/angleCorrectionTest/oldResults --keep ~/data/M112646261LE.IMG ~/data/M112646261RE.IMG
#exit 0

# original file offset = -46, -28


# Move to test directory
#mkdir /home/smcmich1/data/stereoCorrectionTest
#cd /home/smcmich1/data/stereoCorrectionTest



# Copy and convert the first set of cubes
#lronac2isis from=~/data/ARISTARCHU2/M120168714LE.IMG to=M120168714LE.cub
#lronac2isis from=~/data/ARISTARCHU2/M120168714RE.IMG to=M120168714RE.cub





# Convert control point file into usable format
#/home/smcmich1/repo/lronacPipeline/extractQtieControlPoints.py --cnetPath ~/data/ARISTARCHU2/control_pointreg6.net --output controlPoints.csv

# TODO: Apply rotation and offset to second pair (need absolute rot changer?)




# Noproj the corrected data
noproj from=/home/smcmich1/data/stereoCorrectionTest/M120168714LE.final.cub    to=/home/smcmich1/data/stereoCorrectionTest/M120168714LE.final.noproj.cub    match=/home/smcmich1/data/stereoCorrectionTest/M120168714LE.final.cub specs=/home/smcmich1/data/stereoCorrectionTest/noprojInstruments_fullRes.pvl

noproj from=/home/smcmich1/data/stereoCorrectionTest/M120168714RE.final.cub    to=/home/smcmich1/data/stereoCorrectionTest/M120168714RE.final.noproj.cub    match=/home/smcmich1/data/stereoCorrectionTest/M120168714LE.final.cub specs=/home/smcmich1/data/stereoCorrectionTest/noprojInstruments_fullRes.pvl

cp /home/smcmich1/data/stereoCorrectionTest/M120168714LE.final.noproj.cub /home/smcmich1/data/stereoCorrectionTest/M120168714LE.mosaic.cub
handmos from=/home/smcmich1/data/stereoCorrectionTest/M120168714RE.final.noproj.cub mosaic=/home/smcmich1/data/stereoCorrectionTest/M120168714LE.mosaic.cub outsample=0 outline=0 matchbandbin=FALSE priority=ontop

# TODO: Generate DEM, compare to ASU/LOLA




