#!/bin/bash  

# Generate comparison data
#mkdir /home/smcmich1/data/angleCorrectionTest/oldResults
#lronac2mosaic.py -t 4 -o /home/smcmich1/data/angleCorrectionTest/oldResults --keep ~/data/M112646261LE.IMG ~/data/M112646261RE.IMG
#exit 0

# original file offset = -46, -28


# Move to test directory
#mkdir /home/smcmich1/data/stereoCorrectionTest
cd /home/smcmich1/data/stereoCorrectionTest



# Copy and convert the first set of cubes
#lronac2isis from=~/data/ARISTARCHU2/M120168714LE.IMG to=M120168714LE.cub
#lronac2isis from=~/data/ARISTARCHU2/M120168714RE.IMG to=M120168714RE.cub

# Perform pre-processing steps
#lronaccal from=M120168714LE.cub to=M120168714LE.lronaccal.cub
#lronaccal from=M120168714RE.cub to=M120168714RE.lronaccal.cub

#lronacecho from=M120168714LE.lronaccal.cub to=M120168714LE.lronaccal.lronacecho.cub
#lronacecho from=M120168714RE.lronaccal.cub to=M120168714RE.lronaccal.lronacecho.cub

# Apply corrected SPICE data
#/home/smcmich1/repo/lronacPipeline/positionCorrector.py --input ./M120168714LE.lronaccal.lronacecho.cub --output ./M120168714LE.corrected.cub
#/home/smcmich1/repo/lronacPipeline/positionCorrector.py --input ./M120168714RE.lronaccal.lronacecho.cub --output ./M120168714RE.corrected.cub --spk reCorrectedSpk.bsp


# Copy and convert the second set of cubes
#lronac2isis from=~/data/ARISTARCHU2/M120175500LE.IMG to=M120175500LE.cub
#lronac2isis from=~/data/ARISTARCHU2/M120175500RE.IMG to=M120175500RE.cub

# Perform pre-processing steps
#lronaccal from=M120175500LE.cub to=M120175500LE.lronaccal.cub
#lronaccal from=M120175500RE.cub to=M120175500RE.lronaccal.cub

#lronacecho from=M120175500LE.lronaccal.cub to=M120175500LE.lronaccal.lronacecho.cub
#lronacecho from=M120175500RE.lronaccal.cub to=M120175500RE.lronaccal.lronacecho.cub

# Apply corrected SPICE data
#/home/smcmich1/repo/lronacPipeline/positionCorrector.py --input ./M120175500LE.lronaccal.lronacecho.cub --output ./M120175500LE.corrected.cub
#/home/smcmich1/repo/lronacPipeline/positionCorrector.py --input ./M120175500RE.lronaccal.lronacecho.cub --output ./M120175500RE.corrected.cub --spk reCorrectedSpk.bsp

# TODO: Handle kernel-spanning data!


# Convert control point file into usable format
#/home/smcmich1/repo/lronacPipeline/extractQtieControlPoints.py --cnetPath ~/data/ARISTARCHU2/control_pointreg6.net --output controlPoints.csv

# Compute the global rotation and offet between the two LE cubes, recording the GDC stereo points
/home/smcmich1/repo/StereoPipeline/src/asp/Tools/lronacAngleSolver --outputPath ./solvedParams.txt --gdcPointsOutPath stereoGdcPoints.csv --matchingPixelsPath controlPoints.csv ./M120168714LE.corrected.cub ./M120175500LE.corrected.cub --worldTransform --includePosition


# TODO: Call rotation corrector to compute offset (and position?) between two of the images, record 3d points.  Use GLOBAL transform!

# TODO: Apply rotation and offset to second pair (need absolute rot changer?)

# TODO: Use pc-align to compare points to LOLA DEM, compute rotation and offset

# TODO: Apply rotation and offset to LE images

# TODO: Use rotation corrector to rotate RE images



# Noproj the corrected data
#noproj from=M112646261LE.corrected.cub    to=M112646261LE.corrected.noproj.cub    match=M112646261LE.corrected.cub specs=noprojInstruments_fullRes.pvl
#noproj from=M112646261RE.rotCorrected.cub to=M112646261RE.rotCorrected.noproj.cub match=M112646261LE.corrected.cub specs=noprojInstruments_fullRes.pvl 

#cp M112646261LE.corrected.noproj.cub M112646261LE.mosaic.cub
#handmos from=M112646261RE.rotCorrected.noproj.cub mosaic=M112646261LE.mosaic.cub outsample=0 outline=0 matchbandbin=FALSE priority=ontop

# TODO: Generate DEM, compare to ASU/LOLA




