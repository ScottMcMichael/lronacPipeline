#!/bin/bash

# Before running, cd to where you want all this installed and fill in BASE_SYSTEM_DIR

BASE_SYSTEM_DIR=/home/smcmich1/repo/BaseSystem

# Install SuiteSparse
wget http://www.cise.ufl.edu/research/sparse/SuiteSparse/current/SuiteSparse.tar.gz
tar -zxvf SuiteSparse.tar.gz
cd SuitSparse
make
# Do a sort of installation by hand since we can't call make install
mkdir include
find . -name '*.h' -print0  | xargs -0 -I {} cp {} ./include
find . -name '*.hpp' -print0  | xargs -0 -I {} cp {} ./include
mkdir lib
find . -name '*.a' -print0  | xargs -0 -I {} cp {} ./lib 

#TODO: Store the include and lib paths and use below!

# Install CERES
cd ..
git clone https://ceres-solver.googlesource.com/ceres-solver
cd ceres-solver
mkdir build
cd build
cmake -D EIGEN_INCLUDE_DIR_HINTS=$BASE_SYSTEM_DIR/include/eigen3  -D SUITESPARSE_INCLUDE_DIR_HINTS=/home/smcmich1/repo/SuiteSparse/include -D SUITESPARSE_LIBRARY_DIR_HINTS=/home/smcmich1/repo/SuiteSparse/lib -D MINIGLOG=TRUE ..
