#===============================================================================
# Configuration

project(lronacPipeline)
cmake_minimum_required(VERSION 2.8)

# Add path for custom modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}" )

set(PACKAGE_VERSION "0.1")
set(PACKAGE_NAME "LRONAC Pipeline")
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")
set(PACKAGE_BUGREPORT "scott.t.mcmichael@nasa.gov")


#TODO: Remove this, the find functions should handle it!
set(BASESYSTEM_INSTALL_DIR  "/home/smcmich1/repo/BaseSystem")
set(VISIONWORKBENCH_ROOT    "/home/smcmich1/repo/visionworkbench/build")
set(StereoPipeline_ROOT     "/home/smcmich1/repo/StereoPipeline/build")
set(CERES_INSTALL_DIR       "/home/smcmich1/repo/ceres-solver")
set(SUITESPARSE_INSTALL_DIR "/home/smcmich1/repo/SuiteSparse")

# Custom options
option(BUILD_SHARED_LIBS "Produce shared libraries." TRUE)


set(Boost_USE_STATIC_LIBS   OFF)
set(Boost_USE_MULTITHREADED ON )
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})


set(CMAKE_BUILD_TYPE RelWithDebInfo)


#==================================================================================
# Find dependencies

#find_package(Boost 1.50.0
#  COMPONENTS
#  program_options
#  system
#  filesystem
#  thread
#  )

message("--------------------------------------------------------------")
message("Searching for boost...")

set(BOOST_ROOT "${BASESYSTEM_INSTALL_DIR}/include/boost-1_54")
set(BOOST_LIBRARYDIR  "${BASESYSTEM_INSTALL_DIR}/lib")
set(BOOST_INCLUDEDIR "${BASESYSTEM_INSTALL_DIR}/include/boost-1_54/boost")
set(Boost_NO_SYSTEM_PATHS ON) #This does not seem to work
#set(Boost_DEBUG 1) # Uncomment for useful finding boost debug output
set(BOOST_MIN_VERSION "1.54.0")
set(REQUIRED_BOOST_LIBS program_options system filesystem regex date_time thread iostreams)

find_package(Boost ${BOOST_MIN_VERSION} REQUIRED COMPONENTS ${REQUIRED_BOOST_LIBS})

if(Boost_FOUND)
    message("Successfully found Boost version ${Boost_VERSION}!")
    message("---Boost Include Dir = ${Boost_INCLUDE_DIR}")
    #set(Boost_LIBRARY_DIR ${BASESYSTEM_INSTALL_DIR}/lib) # Override this value
    message("---Boost Library Dir = ${Boost_LIBRARY_DIRS}")
    message("---Boost libraries = ${Boost_LIBRARIES}")
    include_directories(${Boost_INCLUDE_DIR})
    link_directories(${Boost_LIBRARY_DIR})
endif(Boost_FOUND)
message("--------------------------------------------------------------")


message("--------------------------------------------------------------")
message("Searching for QT...")
# TODO: Make sure this finds the correct install!
#set (QT_QMAKE_EXECUTABLE  ${BASESYSTEM_INSTALL_DIR}/bin/qmake)
#SET(CMAKE_PREFIX_PATH ${BASESYSTEM_INSTALL_DIR} ${CMAKE_PREFIX_PATH}) 
#set(REQUIRED_QT_LIBS QtCore QtXml QtGui ) #TODO: Figure this out
#find_package(Qt4 4.8.5 REQUIRED )# ${REQUIRED_QT_LIBS})

# Testing: Manual QT settings
set(QT_FOUND True)

set(QT_INCLUDE_DIRS ${BASESYSTEM_INSTALL_DIR}/include/Qt
            ${BASESYSTEM_INSTALL_DIR}/include/Qt3Support
            ${BASESYSTEM_INSTALL_DIR}/include/QtCore
            ${BASESYSTEM_INSTALL_DIR}/include/QtDesigner
            ${BASESYSTEM_INSTALL_DIR}/include/QtGui
            ${BASESYSTEM_INSTALL_DIR}/include/QtHelp
            ${BASESYSTEM_INSTALL_DIR}/include/QtMultimedia
            ${BASESYSTEM_INSTALL_DIR}/include/QtNetwork
            ${BASESYSTEM_INSTALL_DIR}/include/QtSql
            ${BASESYSTEM_INSTALL_DIR}/include/QtSvg
            ${BASESYSTEM_INSTALL_DIR}/include/QtTest
            ${BASESYSTEM_INSTALL_DIR}/include/QtUiTools
            ${BASESYSTEM_INSTALL_DIR}/include/QtXml
            ${BASESYSTEM_INSTALL_DIR}/include/QtXmlPatterns)
set(QT_LIBRARY_DIR  ${BASESYSTEM_INSTALL_DIR}/lib)
set(QT_LIBRARIES  libQt3Support.so
          libQtCLucene.so
          libQtCore.so
          libQtDesignerComponents.so
          libQtDesigner.so
          libQtGui.so
          libQtHelp.so
          libQtMultimedia.so
          libQtNetwork.so
          libQtSql.so
          libQtSvg.so
          libQtTest.so
          libQtXmlPatterns.so
          libQtXml.so )

if(QT_FOUND)
    message("QT include dirs  ")
    foreach(f ${QT_INCLUDE_DIRS}) 
        message("${f}")
        include_directories("${f}")
    endforeach(f) 
    link_directories(${QT_LIBRARY_DIR})
endif(QT_FOUND)
message("--------------------------------------------------------------")

# Get NAIF includes (Needed for reading/writing SPICE data)
set(NAIF_FOUND True)
set(NAIF_INCLUDE_DIRS ${BASESYSTEM_INSTALL_DIR}/include/naif/)



find_package(VisionWorkbench REQUIRED
  COMPONENTS
  core
  math
  image
  camera
  fileio
  cartography
  geometry
  interestPoint
  stereo
  mosaic
  )

message("Found vision workbench include directory: ${VISIONWORKBENCH_INCLUDE_DIRS}")


find_package(StereoPipeline REQUIRED COMPONENTS core spiceIO isisIO)



message("Found stereo pipeline include directory: ${StereoPipeline_INCLUDE_DIR}")


set(PROTOBUF_FOUND)
set(PROTOBUF_INCLUDE_DIR ${BASESYSTEM_INSTALL_DIR}/include   )
set(PROTOBUF_LIBRARY     ${BASESYSTEM_INSTALL_DIR}/lib/libprotobuf.so)
include_directories("${PROTOBUF_INCLUDE_DIR}")
include_directories("${QT_INCLUDE_DIRS}")
include_directories("${NAIF_INCLUDE_DIRS}")


include_directories(${VISIONWORKBENCH_INCLUDE_DIRS})
include_directories(${StereoPipeline_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${PROTOBUF_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

link_directories(${Boost_LIBRARY_DIRS})

set(EIGEN_INCLUDE_DIR ${BASESYSTEM_INSTALL_DIR}/include/eigen3)
include_directories("${EIGEN_INCLUDE_DIR}")

include_directories("${SUITESPARSE_INSTALL_DIR}/include")
link_directories("${SUITESPARSE_INSTALL_DIR}/lib")
set(SUITESPARSE_LIBRARIES  libamd.a   libccolamd.a  libcolamd.a    libklu.a   libspqr.a libbtf.a    libcxsparse.a   libldl.a libsuitesparseconfig.a  libcholmod.a  libcxsparse.a  librbio.a  libumfpack.a libcamd.a)
# The order of the libraries here is important

FIND_PACKAGE(OpenMP REQUIRED)
if (OPENMP_FOUND)

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

#  if (UNIX)
#    # At least on Linux, we need pthreads to be enabled for mutex to
#    # compile.  This may not work on Windows or Android.
#    FIND_PACKAGE(Threads REQUIRED)
#    SET(STATIC_LIBRARY_FLAGS
#      "${STATIC_LIBRARY_FLAGS} ${CMAKE_THREAD_LIBS_INIT}")
#    SET(CMAKE_SHARED_LINKER_FLAGS
#      "${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_THREAD_LIBS_INIT}")
#    ADD_DEFINITIONS(-DCERES_HAVE_PTHREAD)
#    ADD_DEFINITIONS(-DCERES_HAVE_RWLOCK)
#  endif (UNIX)
endif (OPENMP_FOUND)



# Add the Ceres library
set(CERES_LIBRARIES    ${CERES_INSTALL_DIR}/lib/libceres.so ${CERES_INSTALL_DIR}/lib/libminiglog.so)
set(CERES_INCLUDE_DIRS ${CERES_INSTALL_DIR}/include ${CERES_INSTALL_DIR}/include/ceres/internal/miniglog)
foreach(f ${CERES_INCLUDE_DIRS})
    message("${f}")
    include_directories("${f}")
endforeach(f) 




#===================================================================================
# Executables

# TODO: Don't build these, just use the "Tools" repository
add_executable(imagestats imagestats.cc) 
target_link_libraries(imagestats  ${VISIONWORKBENCH_LIBRARIES} ${Boost_LIBRARIES}) 


add_executable(matchBinaryToCsv matchBinaryToCsv.cc) 
target_link_libraries(matchBinaryToCsv  ${VISIONWORKBENCH_LIBRARIES} ${Boost_LIBRARIES}) 

add_executable(pixelPairsFromStereo pixelPairsFromStereo.cc) 
target_link_libraries(pixelPairsFromStereo  ${VISIONWORKBENCH_LIBRARIES}  ${StereoPipeline_LIBRARIES} ${Boost_LIBRARIES}) 


add_executable(lola_compare lola_compare.cc) 
target_link_libraries(lola_compare  ${VISIONWORKBENCH_LIBRARIES} ${Boost_LIBRARIES}) 

add_executable(pixelPairsFromStereo pixelPairsFromStereo.cc) 
target_link_libraries(pixelPairsFromStereo  ${VISIONWORKBENCH_LIBRARIES}  ${StereoPipeline_LIBRARIES}) 




add_executable(stereoIpFind stereoIpFind.cc)
target_link_libraries(stereoIpFind  ${VISIONWORKBENCH_LIBRARIES}  ${StereoPipeline_LIBRARIES}) 

add_executable(lronacAngleDoubleSolver IsisInterfaceLineScanRot.h IsisInterfaceLineScanRot.cc lronacAngleDoubleSolver.cc lronacSolverSupport.h lronacSolverModelDouble.h) 
target_link_libraries(lronacAngleDoubleSolver  ${QT_LIBRARIES} ${VISIONWORKBENCH_LIBRARIES} ${StereoPipeline_LIBRARIES}  ${CERES_LIBRARIES} ${SUITESPARSE_LIBRARIES} ${Boost_LIBRARIES} -lisis3) 

add_executable(spiceEditor IsisInterfaceLineScanRot.h IsisInterfaceLineScanRot.cc SpiceEditor.cc) 
target_link_libraries(spiceEditor  ${VISIONWORKBENCH_LIBRARIES} ${StereoPipeline_LIBRARIES} ${Boost_LIBRARIES} -lisis3 -lcspice) 








