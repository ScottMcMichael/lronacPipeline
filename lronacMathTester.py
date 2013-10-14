#!/usr/bin/env python

import os, glob, optparse, re, shutil, subprocess, sys, string, time

import numpy

def rotZ(angle):
    r = numpy.matrix([[numpy.cos(angle), numpy.sin(angle), 0],[-numpy.sin(angle), numpy.cos(angle), 0],[0, 0, 1]])
    return r

def rotY(angle):
    r = numpy.matrix([[numpy.cos(angle), 0, numpy.sin(angle)],[0, 1, 0],[-numpy.sin(angle), 0, numpy.cos(angle)]])
    return r

def rotX(angle):
    r = numpy.matrix([[1, 0, 0],[0, numpy.cos(angle), numpy.sin(angle)],[0, -numpy.sin(angle), numpy.cos(angle)]])
    return r

def main():

    numpy.set_printoptions(precision=5)
    numpy.set_printoptions(suppress=True)
  
    # Instrument rotation loaded from SPICE (no correction)
    #r_inst = numpy.matrix([[0.725024,  0.368843, -0.581632],[0.1303, -0.902714, -0.410035],[-0.676286, 0.221498, -0.70255]])
    r_inst = numpy.matrix([[0.725015,0.368846,-0.58164],[0.1303,-0.902714,-0.410036],[-0.676295,0.221494,-0.702542]])
    print 'r_inst =  '
    print str(r_inst)
    print ''

    # Spacecraft rotation loaded from SPICE (with LRONAC-R rotation set to 0,0,0 in frame file)
    r_sc = numpy.matrix([[-0.725978, -0.36853, 0.580639 ],[-0.143581,  0.906914,  0.396095 ],[-0.672562,  0.204188, -0.711314 ]])
    print 'r_sc = '
    print str(r_sc)
    print ''

    ## Instrument rotation from original frame file initialized with ZYX function
    #r_zyx = numpy.matrix([[-0.999999, -1.22461e-16, -0.00137881 ],[2.71674e-05, -0.999806, -0.0197035 ],[-0.00137854,  -0.0197035,    0.999805 ]])
    #print 'r_zyx = '
    #print str(r_zyx)
    #print ''

    #TODO: Try this with Y angle negated?
    # Instrument rotation from original frame file initialized with XYZ function
    r_xyz = numpy.matrix ([[-0.999999,2.71674e-05,0.00137854 ], [1.22461e-16,-0.999806,0.0197035 ], [0.00137881,0.0197035,0.999805 ]])
    print 'r_xyz = '
    print str(r_xyz)
    print ''

    print 'Make sure file offset is used correctly - should be zero!'
    print 'numpy.transpose(r_xyz) * r_sc - (r_inst)'
    print numpy.transpose(r_xyz) * r_sc - (r_inst)
    print ''


    ## Offset matrix computed by solution (X,Y only)
    #r_offsetXY = numpy.matrix([[1, 1.71783e-07, -0.00037477 ], [0, 1, 0.000458369 ], [0.00037477, -0.000458369,  1 ]])
    #print 'r_offsetXY = '
    #print str(r_offsetXY)
    #print ''

    # Offset matrix computed by solution (X,Y,Z)
    #r_offsetXYZ = numpy.matrix([[0.999976, -0.00697023, -0.000199261 ], [0.00697032, 0.999976, 0.00045699 ], [0.000196071, -0.000458368, 1 ]])
    r_offsetXYZ = numpy.matrix([[0.999997,0.00243865,5.3149e-05 ], [-0.00243868,0.999997,0.000458199], [-5.20315e-05,-0.000458327,1]])
    print 'r_offsetXYZ = '
    print str(r_offsetXYZ)
    print '\n'


    ## Check the two methods of arriving at the same result
    #print 'XY solution error: '
    #print numpy.transpose(r_xyz * numpy.transpose(r_offsetXY)) * r_sc - (r_offsetXY*r_inst)
    #print ''
    
    print 'XYZ solution error (should be zero): '
    print numpy.transpose(r_xyz * numpy.transpose(r_offsetXYZ)) * r_sc - (r_offsetXYZ*r_inst)
    print ''

    # Read in from modified frame file with xyz function
    r_xyz_modified = numpy.matrix([[-0.999996,0.0024728,0.00143043], [-0.00244346,-0.999794,0.0201618], [0.00148,0.0201582,0.999796]])

    #r_xzy_modified = numpy.matrix([[-0.999996,0.00241313,-0.00152895], [-0.00244346,-0.999794,0.0201582], [-0.00147999,0.0201618,0.999796]])


    print 'r_offsetXYZ*r_inst'
    print r_offsetXYZ*r_inst
    print ''

    #print 'print r_xyz*numpy.transpose(r_offsetXYZ)'
    #print r_xyz*numpy.transpose(r_offsetXYZ)
    #print ''

    print 'transpose(r_xyz_modified) * r_sc'
    print numpy.transpose(r_xyz_modified) * r_sc
    print ''

    #print 'transpose(r_xzy_modified) * r_sc'
    #print numpy.transpose(r_xzy_modified) * r_sc
    #print ''

    print 'numpy.transpose(r_xyz * numpy.transpose(r_offsetXYZ)) * r_sc'
    print numpy.transpose(r_xyz * numpy.transpose(r_offsetXYZ)) * r_sc
    print ''

    print 'R_inst (getting) = Matrix3x3((0.724626,0.371074,-0.580708)(0.131764,-0.901709,-0.411776)(-0.676428,0.221867,-0.702296))'
    print 'R_inst (correct) = Matrix3x3((0.725292,0.366656,-0.582678)(0.128221,-0.903509,-0.408938)(-0.676395,0.221888,-0.702322))'
    print ''

    # This method lets us just stick a rotation matrix in the offset file
    # We want it to be the eqivalent of using this offset matrix on the original instrument base: Matrix3x3((0.999997,0.00243865,5.3149e-05)(-0.00243868,0.999997,0.000458199)(-5.20315e-05,-0.000458327,1))
    r_put_in_file = numpy.matrix([[-0.999996, 0.00246648, 0.00143056], [-0.00243713, -0.999794, 0.0201617], [0.00147999, 0.0201582, 0.999796]])
    print 'r_put_in_file = '
    print str(r_put_in_file)
    print ''

    r_put_in_fileT = numpy.matrix([[-0.999996,-0.00243713,0.00147999], [0.00246648,-0.999794,0.0201582], [0.00143056,0.0201617,0.999796]])
    print 'r_put_in_fileT = '
    print str(r_put_in_fileT)
    print ''

    print 'numpy.transpose(r_put_in_file) * r_sc'
    print numpy.transpose(r_put_in_file) * r_sc
    print ''

    print 'numpy.transpose(r_put_in_fileT) * r_sc'
    print numpy.transpose(r_put_in_fileT) * r_sc
    print ''

    return 0


    print 'SC rotation by 30 deg x, y, z'
    deg2rad = 3.14159265359/180
    print numpy.transpose(rotX(30*deg2rad)) * r_sc
    print numpy.transpose(rotY(-30*deg2rad)) * r_sc
    print numpy.transpose(rotZ(30*deg2rad)) * r_sc

#R_inst   = Matrix3x3((-0.725978,-0.36853,0.580639)(0.211937,0.683317,0.698685)(-0.654247,0.630289,-0.417968))
#R_inst   = Matrix3x3((-0.969173,-0.209879,0.129052)(-0.143672,0.90694,0.396002)(-0.200155,0.365254,-0.909136))
#R_inst   = Matrix3x3((-0.556971,-0.772599,0.304753) (-0.487361,0.601155,0.633318) (-0.672505,0.204215,-0.711361))

    print 'rot 10 0 20'
    print numpy.transpose( rotX(10*deg2rad)*rotZ(20*deg2rad)) * r_sc

#R_inst   = Matrix3x3((-0.675069,-0.639181,0.368415)(-0.272182,0.679923,0.680897)(-0.685711,0.359376,-0.632969))


#Y rotation appears to be negated between vision workbench definition and ISIS treatment!
#Does this need to be handled in the C code?
#Try euler back-calc function with Y rotation matrix sin signs flipped 

    print 'rot 10 20 0'
    print numpy.transpose( rotX(10*deg2rad)*rotY(20*deg2rad)) * r_sc


    print numpy.transpose( rotX(10*deg2rad)*rotY(-20*deg2rad)) * r_sc

#R_inst   = Matrix3x3((-0.921958,-0.21978,0.31889)(-0.026216,0.856913,0.514794)(-0.386402,0.466258,-0.795799))

    print 'rot 10 20 30'
    print numpy.transpose( rotX(10*deg2rad)*rotY(20*deg2rad)*rotZ(30*deg2rad) ) * r_sc
#    print numpy.transpose( rotZ(30*deg2rad)*rotY(20*deg2rad)*rotX(10*deg2rad) ) * r_sc
#    print numpy.transpose( rotX(30*deg2rad)*rotY(20*deg2rad)*rotZ(10*deg2rad) ) * r_sc
#    print numpy.transpose( rotZ(10*deg2rad)*rotY(20*deg2rad)*rotX(30*deg2rad) ) * r_sc

    print numpy.transpose( rotX(10*deg2rad)*rotY(-20*deg2rad)*rotZ(30*deg2rad) ) * r_sc

#R_inst   = Matrix3x3((-0.786303,-0.61764,0.0157479)(-0.484828,0.632622,0.603929)(-0.382973,0.467236,-0.796882))


    print '\nrot -0.026260, 0.002981, -0.139726'
    # 1.15506, 0.0847975, -179.86 in file
    print numpy.transpose(r_xyz * rotX(1.15506*deg2rad)*rotY(0.0847975*deg2rad)*rotZ(-179.86*deg2rad) ) * r_sc
    print numpy.transpose(r_xyz * rotX(1.15506*deg2rad)*rotY(-0.0847975*deg2rad)*rotZ(-179.86*deg2rad) ) * r_sc
    print numpy.transpose(r_xyz * rotX(1.15506*deg2rad)*rotY(0*deg2rad)*rotZ(-179.86*deg2rad) ) * r_sc
    #print numpy.transpose(r_xyz * rotX(-0.026260*deg2rad)*rotY(-0.002981*deg2rad)*rotZ(-0.139726*deg2rad) ) * r_sc

# Current: R_inst   = Matrix3x3((0.726691,0.370394,-0.578557)(0.131771,-0.901712,-0.411767)(-0.674208,0.22299,-0.704073))

    # This output matches the computed combined pose
    #print numpy.transpose(r_xyz * rotX(-0.026260*deg2rad)*rotY(0.002981*deg2rad)*rotZ(-0.139726*deg2rad) ) * r_sc



if __name__ == "__main__":
    sys.exit(main())