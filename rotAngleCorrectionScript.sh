#!/bin/bash  

# Generate comparison data
#mkdir /home/smcmich1/data/angleCorrectionTest/oldResults
#lronac2mosaic.py -t 2 -o /home/smcmich1/data/angleCorrectionTest/oldResults --keep ~/data/LinneCrater/M139836046LE.IMG ~/data/LinneCrater/M139836046RE.IMG
#exit 0

# original file offset = -46, -28


# Move to test directory
mkdir /home/smcmich1/data/angleCorrectionTest
cd /home/smcmich1/data/angleCorrectionTest


# Copy and convert the two test cubes
#lronac2isis from=~/data/LinneCrater/M139836046LE.IMG to=M139836046LE.cub
#lronac2isis from=~/data/LinneCrater/M139836046RE.IMG to=M139836046RE.cub

# Perform pre-processing steps
#lronaccal from=M139836046LE.cub to=M139836046LE.lronaccal.cub
#lronaccal from=M139836046RE.cub to=M139836046RE.lronaccal.cub

#lronacecho from=M139836046LE.lronaccal.cub to=M139836046LE.lronaccal.lronacecho.cub
#lronacecho from=M139836046RE.lronaccal.cub to=M139836046RE.lronaccal.lronacecho.cub


# Apply corrected SPICE data
/home/smcmich1/repot/lronacPipeline/positionCorrector.py --input M139836046LE.lronaccal.lronacecho.cub --output M139836046LE.corrected.cub
/home/smcmich1/repot/lronacPipeline/positionCorrector.py --input M139836046RE.lronaccal.lronacecho.cub --output M139836046RE.corrected.cub

# TODO: Handle kernel-spanning data!

# Noproj the corrected data
noproj from=M139836046LE.corrected.cub to=M139836046LE.corrected.noproj.cub match=M139836046LE.corrected.cub specs=~/repot/lronacPipeline/test.pvl
noproj from=M139836046RE.corrected.cub to=M139836046RE.corrected.noproj.cub match=M139836046LE.corrected.cub specs=~/repot/lronacPipeline/test.pvl

# TODO: evaluate offset

# TODO: Least squares solver C++ program to evaluate angle between cameras
# - Stereo::StereoModel
# - asp::IsisIo::IsisCameraModel::point_to_pixel
# - Math::LeastSquaresBase


# TODO: Apply the rotation to the cameras

# TODO: Generate DEM, compare to ASU/LOLA




