#
# source /nusoft/app/externals/setup
# setup cmake
# setup genie  v2_8_0a   -q debug:e4
# setup root   v5_34_05  -q debug:e4
#
# cd /nova/app/users/$USER
#
# mkdir dk2nu-build dk2n-install
# cd    dk2nu-build
#
# cmake [ -DWITH_GENIE=ON ] \
#       [ -DGENIE_ONLY=OFF -DDK2NUDATA_DIR=<path-to-dk2nudata-install> ] \
#       [ -DWITH_TBB=OFF ] \
#     -DCMAKE_INSTALL_PREFIX:PATH=<path-to-install>/dk2nu-install \
#     <path-to-source>/dk2nu
#
# Required environment variables:
#    ROOTSYS    e.g. /nusoft/app/externals/root/v5_34_05/Linux64bit+2.6-2.5-e2-debug
# Optional
#    GENIE
#      if GENIE, also requires  $LIBXML2_INC and $LOG4CPP_INC
#
#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(dk2nu)
set(CMAKE_BUILD_TYPE Debug)
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

#----------------------------------------------------------------------------
# By default attempt to build GENIE flux driver
# You can set WITH_GENIE to OFF via the command line or ccmake/cmake-gui
# to build only the base class library
#
option(WITH_GENIE "Build GENIE flux driver" ON)
option(GENIE_ONLY "Build only GENIE flux driver" OFF)
option(DK2NUDATA_DIR "path the dk2nu data half of package if GENIE_ONLY" "")
option(WITH_TBB "Build in linkage to TBB library" ON)
option(COPY_AUX "install etc, convert, snippets subdirectories" ON)

execute_process(COMMAND root-config --version OUTPUT_VARIABLE ROOT_VERSION)
string(STRIP ${ROOT_VERSION} ROOT_VERSION)
string(SUBSTRING ${ROOT_VERSION} 0 1 ROOT_MAJOR_VERSION)
MESSAGE("--DK2NU-- ROOT_VERSION=${ROOT_VERSION} major=${ROOT_MAJOR_VERSION}")
# for now this still has to work with both 5 & 6

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake
                      $ENV{ROOTSYS}/cmake/modules
                      $ENV{ROOTSYS}/etc/cmake
                      ${CMAKE_MODULE_PATH} )

#MESSAGE("CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")
#MESSAGE(" ")

find_package(ROOT REQUIRED)

include(Darwin)

if( ${ROOT_MAJOR_VERSION} LESS 6 )
  # older ROOT wasn't so helpful, mimic what ROOT 6 does

  MESSAGE("--DK2NU-- start ROOT 5 config")
  include(CMakeParseArguments)
  ##  $ENV{ROOTSYS}/etc/cmake/ROOTUseFile.cmake)

  execute_process(COMMAND root-config --auxcflags
                  OUTPUT_VARIABLE ROOT_CXX_FLAGS)
  ### has a troublesome trailing \n
  string(STRIP ${ROOT_CXX_FLAGS} ROOT_CXX_FLAGS)
  MESSAGE("--DK2NU-- ROOT_CXX_FLAGS ${ROOT_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS})

  include(CMakeParseArguments)

  find_library(ROOT_Physics_LIBRARY NAMES Physics PATHS
    $ENV{ROOTSYS}/lib NO_DEFAULT_PATH )
  find_library(ROOT_Matrix_LIBRARY NAMES Matrix PATHS
    $ENV{ROOTSYS}/lib NO_DEFAULT_PATH )

else()
  MESSAGE("--DK2NU-- start ROOT 6 config")
  # don't know who defines the ${ROOT_DIR} used here
  # seems to be $ENV{ROOTSYS}/cmake
  include(${ROOT_DIR}/ROOTUseFile.cmake)
endif()
MESSAGE("--DK2NU-- finish ROOT X config")

if( WITH_TBB )
  MESSAGE("--DK2NU-- enabled WITH_TBB")
  find_library(TBB_LIBRARY NAMES tbb PATHS ENV TBB_LIB NO_DEFAULT_PATH )
  # ROOT_Imt_LIBRARY = implicit multi-threading
  set(ADD_ROOT_Imt_LIBRARY ${ROOT_Imt_LIBRARY})
else()
  MESSAGE("--DK2NU-- WITH_TBB disabled")
  set(TBB_LIBRARY "")
  set(ADD_ROOT_Imt_LIBRARY "")
endif()

#----------------------------------------------------------------------------

if( ${ROOT_MAJOR_VERSION} LESS 6 )
  #### / ) # add abs path for find_file() in FindROOT.cmake
  set(ROOT5_NEEDS_FULL_PATH / )
endif()


MESSAGE("====> PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}")
set(INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/..
                        ${PROJECT_SOURCE_DIR}
                        ${PROJECT_SOURCE_DIR}/dk2nu
                        $${ROOT_INCLUDE_DIR}
                        ${ROOT5_NEEDS_FULL_PATH} )

  set(GENIE $ENV{GENIE})
if( GENIE_ONLY )
  set(COPY_AUX OFF)
  add_subdirectory(genie)
  MESSAGE("--DK2NU-- Disabled dk2nu tree (dk2nudata) building")
else()
  add_subdirectory(tree)
  if ( NOT WITH_GENIE )
    MESSAGE("--DK2NU-- Disabled GENIE flux driver building")
  else()
    add_subdirectory(genie)
    MESSAGE("--DK2NU-- Both tree & GENIE flux driver building enabled")
  endif()
endif()

#----------------------------------------------------------------------------

# Install support files and scripts
if(COPY_AUX)
#  install(FILES README DESTINATION doc RENAME dk2nu.README)
  install(FILES       etc/locations.txt
                      etc/NuMIlocations.txt
                      etc/LBNElocations.txt
          DESTINATION etc)
  install(FILES       scripts/load_dk2nu.C
          DESTINATION scripts)
  install(FILES       scripts/convert/common_convert.C
          DESTINATION scripts/convert)
  install(FILES       scripts/convert/aux/mkgclasses3.sh
          DESTINATION scripts/convert/aux)
  install(FILES       scripts/convert/g3numi/g3numi.C
                      scripts/convert/g3numi/g3numi.h
          DESTINATION scripts/convert/g3numi)
  install(FILES       scripts/convert/g4numi/g4numi.C
                      scripts/convert/g4numi/g4numi.h
          DESTINATION scripts/convert/g4numi)
  install(FILES       scripts/convert/g4minerva/g4minerva.C
                      scripts/convert/g4minerva/g4minerva.h
                      scripts/convert/g4minerva/convert_g4minerva.C
          DESTINATION scripts/convert/g4minerva)
  install(FILES       scripts/convert/flugg/flugg.C
                      scripts/convert/flugg/flugg.h
                      scripts/convert/flugg/convert_flugg.C
          DESTINATION scripts/convert/flugg)
  install(FILES       scripts/convert/g4lbne/g4lbne.C
                      scripts/convert/g4lbne/g4lbne.h
                      scripts/convert/g4lbne/convert_g4lbne.C
          DESTINATION scripts/convert/g4lbne)
  install(FILES       scripts/examples/test_read_locations.C
                      scripts/examples/test_fill_dk2nu.C
                      scripts/examples/test_read_dk2nu.C
          DESTINATION scripts/examples)
endif()
