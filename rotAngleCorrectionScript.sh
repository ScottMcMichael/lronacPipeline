#!/bin/bash  

# Generate comparison data
#mkdir /home/smcmich1/data/angleCorrectionTest/oldResults
#lronac2mosaic.py -t 4 -o /home/smcmich1/data/angleCorrectionTest/oldResults --keep ~/data/M139829261LE.IMG ~/data/M139829261RE.IMG
#exit 0

# original file offset = -46, -28


# Move to test directory
#mkdir /home/smcmich1/data/angleCorrectionTest
cd /home/smcmich1/data/angleCorrectionTest


# Copy and convert the two test cubes
#lronac2isis from=~/data/M139829261LE.IMG to=M139829261LE.cub
#lronac2isis from=~/data/M139829261RE.IMG to=M139829261RE.cub

# Perform pre-processing steps
#lronaccal from=M139829261LE.cub to=M139829261LE.lronaccal.cub
#lronaccal from=M139829261RE.cub to=M139829261RE.lronaccal.cub

#lronacecho from=M139829261LE.lronaccal.cub to=M139829261LE.lronaccal.lronacecho.cub
#lronacecho from=M139829261RE.lronaccal.cub to=M139829261RE.lronaccal.lronacecho.cub

# Apply corrected SPICE data
#/home/smcmich1/repo/lronacPipeline/positionCorrector.py --input M139829261LE.lronaccal.lronacecho.cub --output M139829261LE.corrected.cub
#/home/smcmich1/repo/lronacPipeline/positionCorrector.py --input M139829261RE.lronaccal.lronacecho.cub --output M139829261RE.corrected.cub --spk reCorrectedSpk.bsp

# TODO: Handle kernel-spanning data!

# Now compute the corrected rotation
#/home/smcmich1/repo/lronacPipeline/rotationCorrector.py --left M139829261LE.corrected.cub --right M139829261RE.corrected.cub --output M139829261RE.rotCorrected.cub --spk reCorrectedSpk.bsp --keep

# Noproj the corrected data
#noproj from=M139829261LE.corrected.cub    to=M139829261LE.corrected.noproj.cub    match=M139829261LE.corrected.cub specs=noprojInstruments_fullRes.pvl
#noproj from=M139829261RE.rotCorrected.cub to=M139829261RE.rotCorrected.noproj.cub match=M139829261LE.corrected.cub specs=noprojInstruments_fullRes.pvl 

cp M139829261LE.corrected.noproj.cub M139829261LE.mosaic.cub
handmos from=M139829261RE.rotCorrected.noproj.cub mosaic=M139829261LE.mosaic.cub outsample=0 outline=0 matchbandbin=FALSE priority=ontop

# TODO: Generate DEM, compare to ASU/LOLA




