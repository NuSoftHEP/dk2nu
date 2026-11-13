#
# source /nusoft/app/externals/setup
# setup cmake
# setup genie  <version, R-3 or later>
# setup root   <version, 6 or later>
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
#    ROOTSYS    the ROOT (6 or later) installation
# Optional
#    GENIE
#      if GENIE, also requires  $LIBXML2_INC and $LOG4CPP_INC
#
#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(dk2nu)
# default to an optimised build with debug symbols (-O2 -g);
# -DCMAKE_BUILD_TYPE=<type> on the command line takes precedence
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()
MESSAGE("--DK2NU-- CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
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
option(DK2NU_PYTHON_ENABLED "Enable the building of python bindings" OFF)

execute_process(COMMAND root-config --version OUTPUT_VARIABLE ROOT_VERSION)
string(STRIP ${ROOT_VERSION} ROOT_VERSION)
string(SUBSTRING ${ROOT_VERSION} 0 1 ROOT_MAJOR_VERSION)
MESSAGE("--DK2NU-- ROOT_VERSION=${ROOT_VERSION} major=${ROOT_MAJOR_VERSION}")
if( ${ROOT_MAJOR_VERSION} LESS 6 )
  MESSAGE(FATAL_ERROR "--DK2NU-- ROOT 6 or later is required (found ${ROOT_VERSION})")
endif()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake
                      $ENV{ROOTSYS}/cmake/modules
                      $ENV{ROOTSYS}/etc/cmake
                      ${CMAKE_MODULE_PATH} )

#MESSAGE("CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")
#MESSAGE(" ")

find_package(ROOT REQUIRED)

include(Darwin)

# ${ROOT_DIR} is set by find_package(ROOT)
include(${ROOT_DIR}/ROOTUseFile.cmake)

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

MESSAGE("====> PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}")
set(INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/..
                        ${PROJECT_SOURCE_DIR}
                        ${PROJECT_SOURCE_DIR}/dk2nu
                        ${ROOT_INCLUDE_DIR} )

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

if(DK2NU_PYTHON_ENABLED)
  add_subdirectory(python)
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
