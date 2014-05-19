#
# source /nusoft/app/externals/setup
# setup cmake
# setup genie  v2_8_0a   -q debug:e4
# setup root   v5_34_05  -q debug:e4
#
# cd /nova/app/users/$USER
#
# mkdir dk2nu-build
# cd    dk2nu-build
#
# cmake [ -DWITH_GENIE=OFF ] ../dk2nu # or whatever path
#
# Required environment variables:
#    ROOTSYS    e.g. /nusoft/app/externals/root/v5_34_05/Linux64bit+2.6-2.5-e2-debug
# Optional
#    GENIE
#      if GENIE requires  $LIBXML2_INC and $LOG4CPP_INC
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
option(COPY_AUX "install etc, convert, snippets subdirectories" ON)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake
                      $ENV{ROOTSYS}/cmake/modules
                      $ENV{ROOTSYS}/../source/root/cmake/modules
                      $ENV{ROOTSYS}/etc/cmake 
                      ${CMAKE_MODULE_PATH})

#MESSAGE("CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")
#MESSAGE(" ")

find_package(ROOT REQUIRED)

set(INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/..
                        ${PROJECT_SOURCE_DIR}
                        ${PROJECT_SOURCE_DIR}/tree
                        ${PROJECT_SOURCE_DIR}/genie
                        ${ROOT_INCLUDE_DIR}
                        / ) # add abs path for find_file() in FindROOT.cmake

set(LINK_DIRECTORIES ${ROOT_LIBRARY_DIR} ${CMAKE_CURRENT_BINARY_DIR}/lib )
#MESSAGE("CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}")

if(WITH_GENIE)
  set(GENIE $ENV{GENIE})
  MESSAGE("-- WITH_GENIE was true")
  set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} 
                          ${PROJECT_SOURCE_DIR}/genie
                          $ENV{LIBXML2_INC}
                          $ENV{LOG4CPP_INC}
                          ${GENIE}/src)
  set(LINK_DIRECTORIES ${LINK_DIRECTORIES} ${GENIE}/lib)
  #MESSAGE("--DK2NU- INCLUDE_DIRECTORIES=${INCLUDE_DIRECTORIES}")
endif()

include_directories( ${INCLUDE_DIRECTORIES} )
link_directories( ${LINK_DIRECTORIES} )

#----------------------------------------------------------------------------

include(Darwin)

#include(AltGenDict2)
#include(RootNewMacros)

# $ROOTSYS/etc/cmake/FindROOT.cmake
#    needs CMAKE_PARSE_ARGUMENTS defined to do ROOT_GENERATE_DICTIONARY
# where should one find _that_ 

include(CMakeParseArguments)

#MESSAGE("--DK2NU- PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}")
#MESSAGE("--DK2NU- ROOT_LIBRARIES=${ROOT_LIBRARIES}")
#MESSAGE("--DK2NU-")

#----------------------------------------------------------------------------
#
# libdk2nuTree.so
#
#MESSAGE("--DK2NU- dk2nuTree section begin")
file(GLOB dk2nuTree_SRCS ${PROJECT_SOURCE_DIR}/tree/*.cc 
                         ${PROJECT_SOURCE_DIR}/tree/*.cxx)
file(GLOB dk2nuTree_HDRS ${PROJECT_SOURCE_DIR}/tree/*.h)
list(REMOVE_ITEM dk2nuTree_HDRS ${PROJECT_SOURCE_DIR}/tree/LinkDef.h)

#MESSAGE("dk2nuTree_SRCS=${dk2nuTree_SRCS}")
#MESSAGE(" ")
#MESSAGE("dk2nuTree_HDRS=${dk2nuTree_HDRS}")
#MESSAGE(" ")

set(dk2nuTree_LINKDEF ${PROJECT_SOURCE_DIR}/tree/LinkDef.h)
# do not include ".cxx" here!!!!
set(dk2nuTree_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/dk2nuTreeDict)

#MESSAGE("--DK2NU- about to ROOT_GENERATE_DICTIONARY ${dk2nuTree_DICTIONARY}")
#MESSAGE("--rwh- ROOT_GENERATE_DICTIONARY ${dk2nuTree_DICTIONARY} ${dk2nuTree_HDRS} ${dk2nuTree_LINKDEF} ")
ROOT_GENERATE_DICTIONARY("${dk2nuTree_DICTIONARY}" "${dk2nuTree_HDRS}" LINKDEF "${dk2nuTree_LINKDEF}" )

# add dictionary to list of source files
SET(dk2nuTree_SRCS ${dk2nuTree_SRCS} ${dk2nuTree_DICTIONARY})

add_library(dk2nuTree SHARED ${dk2nuTree_SRCS})
target_link_libraries(dk2nuTree ${ROOT_LIBRARIES} -lPhysics -lMatrix )

#MESSAGE("--DK2NU- dk2nuTree section done")
#MESSAGE(" ")

#----------------------------------------------------------------------------
#
# libdk2nuGenie.so
#
if(WITH_GENIE)

#MESSAGE("--DK2NU- dk2nuTree section begin")
file(GLOB dk2nuGenie_SRCS ${PROJECT_SOURCE_DIR}/genie/*.cc 
                          ${PROJECT_SOURCE_DIR}/genie/*.cxx)
file(GLOB dk2nuGenie_HDRS ${PROJECT_SOURCE_DIR}/genie/*.h)
list(REMOVE_ITEM dk2nuGenie_HDRS ${PROJECT_SOURCE_DIR}/genie/LinkDef.h)

#MESSAGE("dk2nuGenie_SRCS=${dk2nuGenie_SRCS}")
#MESSAGE(" ")
#MESSAGE("dk2nuGenie_HDRS=${dk2nuGenie_HDRS}")
#MESSAGE(" ")

set(dk2nuGenie_LINKDEF ${PROJECT_SOURCE_DIR}/genie/LinkDef.h)
# do not include ".cxx" here!!!!
set(dk2nuGenie_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/dk2nuGenieDict)

#MESSAGE("--DK2NU- about to ROOT_GENERATE_DICTIONARY ${dk2nuGenie_DICTIONARY}")
ROOT_GENERATE_DICTIONARY("${dk2nuGenie_DICTIONARY}" "${dk2nuGenie_HDRS}" LINKDEF "${dk2nuGenie_LINKDEF}" )

# add dictionary to list of source files
SET(dk2nuGenie_SRCS ${dk2nuGenie_SRCS} ${dk2nuGenie_DICTIONARY})

# don't use -l for dk2nuTree if we want dk2nuGenie to depend on dk2nuTree
# before dk2nuTree is built
add_library(dk2nuGenie SHARED ${dk2nuGenie_SRCS})
target_link_libraries(dk2nuGenie ${ROOT_LIBRARIES} -lPhysics -lMatrix dk2nuTree )

#set_property(TARGET dk2nuGenie PROPERTY DEPENDS dk2nuTree)
#add_custom_command(OUTPUT dk2nuGenie COMMAND echo hey there DEPENDS dk2nuTree)

#MESSAGE("--DK2NU- dk2nuGenie section done")

endif()

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
##set(EXAMPLEG4NUMI_SCRIPTS
## me000z200i3.mac
##  )
##
##foreach(_script ${EXAMPLEG4NUMI_SCRIPTS})
##  configure_file(
##    ${PROJECT_SOURCE_DIR}/${_script}
##    ${PROJECT_BINARY_DIR}/${_script}
##    COPYONLY
##    )
##endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
##add_custom_target(G4NUMI DEPENDS g4numi)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS dk2nuTree DESTINATION lib)
if(WITH_GENIE)
  install(TARGETS dk2nuGenie DESTINATION lib)
endif()
#--------------------
# Install the headers
install(FILES ${dk2nuTree_HDRS} DESTINATION include/dk2nu/tree)
if(WITH_GENIE)
  install(FILES ${dk2nuGenie_HDRS} DESTINATION include/dk2nu/genie)
endif()
# Install support files and scripts
if(COPY_AUX)
#  install(FILES README DESTINATION doc RENAME dk2nu.README)
  install(FILES etc/locations.txt etc/LBNElocations.txt DESTINATION etc)
  install(FILES convert/common_convert.C   DESTINATION convert)
  install(FILES convert/aux/mkgclasses3.sh DESTINATION convert/aux)
  install(FILES       convert/g3numi/g3numi.C 
                      convert/g3numi/g3numi.h 
          DESTINATION convert/g3numi)
  install(FILES       convert/g4numi/g4numi.C 
                      convert/g4numi/g4numi.h 
          DESTINATION convert/g4numi)
  install(FILES       convert/g4lbne/g4lbne.C 
                      convert/g4lbne/g4lbne.h 
                      convert/g4lbne/convert_g4lbne.C
          DESTINATION convert/g4lbne)
  install(FILES       convert/g4minerva/g4minerva.C 
                      convert/g4minerva/g4minerva.h 
                      convert/g4minerva/convert_g4minerva.C
          DESTINATION convert/g4minerva)
  install(FILES       convert/flugg/flugg.C 
                      convert/flugg/flugg.h 
                      convert/flugg/convert_flugg.C
          DESTINATION convert/flugg)
  install(FILES       snippets/load_dk2nu.C
                      snippets/test_read_locations.C
                      snippets/test_fill_dk2nu.C
                      snippets/test_read_dk2nu.C
          DESTINATION snippets)
endif()
