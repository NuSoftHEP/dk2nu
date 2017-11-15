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
# cmake [ -DWITH_GENIE=ON ] [ -DWITH_TBB=OFF ]
#     -DCMAKE_INSTALL_PREFIX=../dk2nu-install ../dk2nu # or whatever path
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

set(INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR}/..
                        ${PROJECT_SOURCE_DIR}
                        ${PROJECT_SOURCE_DIR}/tree
                        ${PROJECT_SOURCE_DIR}/dk2nu/tree
                        ${ROOT_INCLUDE_DIR}
                        ${ROOT5_NEEDS_FULL_PATH} )

MESSAGE("--DK2NU-- PROJECT_SOURCE_DIR ${PROJECT_SOURCE_DIR}")
MESSAGE("--DK2NU-- INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES}")


set(LINK_DIRECTORIES ${ROOT_LIBRARY_DIR} ${CMAKE_CURRENT_BINARY_DIR}/lib )
#MESSAGE("CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}")

if(WITH_GENIE)
  set(GENIE $ENV{GENIE})
  MESSAGE("--DK2NU-- WITH_GENIE was true")
  set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES}
                          ${PROJECT_SOURCE_DIR}/genie
                          ${PROJECT_SOURCE_DIR}/dk2nu/genie
                          $ENV{LIBXML2_INC}
                          $ENV{LOG4CPP_INC}
                          ${GENIE}/src)
  set(LINK_DIRECTORIES ${LINK_DIRECTORIES} ${GENIE}/lib)
  #MESSAGE("--DK2NU- INCLUDE_DIRECTORIES=${INCLUDE_DIRECTORIES}")
endif()

include_directories( ${INCLUDE_DIRECTORIES} )
link_directories( ${LINK_DIRECTORIES} )

#----------------------------------------------------------------------------

# ROOT dictionary options;
# GENIE driver dictionary needed "-p" w/ ROOT 5


#MESSAGE("--DK2NU- PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}")
#MESSAGE("--DK2NU- ROOT_LIBRARIES=${ROOT_LIBRARIES}")
#MESSAGE("--DK2NU- CMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}")
#MESSAGE("--DK2NU-")

#----------------------------------------------------------------------------
#
# libdk2nuTree.so
#
#MESSAGE("--DK2NU- dk2nuTree section begin")
file(GLOB dk2nuTree_SRCS ${PROJECT_SOURCE_DIR}/tree/*.cc
                         ${PROJECT_SOURCE_DIR}/tree/*.cxx)

#file(GLOB dk2nuTree_HDRS ${PROJECT_SOURCE_DIR}/tree/*.h)
#list(REMOVE_ITEM dk2nuTree_HDRS ${PROJECT_SOURCE_DIR}/tree/LinkDef.h)
set(dk2nuTree_HDRS dk2nu/tree/dk2nu.h
                   dk2nu/tree/dkmeta.h
                   dk2nu/tree/NuChoice.h
                   dk2nu/tree/dflt.h
                   dk2nu/tree/calcLocationWeights.h
                   dk2nu/tree/readWeightLocations.h )

#MESSAGE("dk2nuTree_SRCS=${dk2nuTree_SRCS}")
#MESSAGE(" ")
#MESSAGE("dk2nuTree_HDRS=${dk2nuTree_HDRS}")
#MESSAGE(" ")

set(dk2nuTree_LINKDEF ${PROJECT_SOURCE_DIR}/tree/LinkDef.h)
# do not include ".cxx" here!!!!
set(dk2nuTree_DICTIONARY dk2nuTreeDict)

#MESSAGE("--DK2NU- about to ROOT_GENERATE_DICTIONARY ${dk2nuTree_DICTIONARY}")
#MESSAGE("--rwh- ROOT_GENERATE_DICTIONARY ${dk2nuTree_DICTIONARY} ${dk2nuTree_HDRS} ${dk2nuTree_LINKDEF} ")
if( ${ROOT_MAJOR_VERSION} LESS 6 )
  set(ROOT_DICT_OPT_TREE -p)
  ROOT_GENERATE_DICTIONARY("${dk2nuTree_DICTIONARY}" "${dk2nuTree_HDRS}"
    LINKDEF "${dk2nuTree_LINKDEF}" OPTIONS ${ROOT_DICT_OPT_TREE} )
else()
  set(ROOT_INCLUDE_PATH $ENV{ROOT_INCLUDE_PATH}:${PROJECT_SOURCE_DIR}:${PROJECT_SOURCE_DIR}/tree)
  # these are defined for the install stage
  set(TREEDICTFILES ${CMAKE_BINARY_DIR}/libdk2nuTree_rdict.pcm
    ${CMAKE_BINARY_DIR}/libdk2nuTree.rootmap )
  set(ROOT_DICT_OPT_TREE -I.. -I../tree
    -I../include -I../include/dk2nu
    -Idk2nu -I../include/tree -I../include/dk2nu/tree -Idk2nu/tree
    -p ) ### -inlineInputHeader )
  # new MODULE option for this funciton to use w/ CMake
  ROOT_GENERATE_DICTIONARY("${dk2nuTree_DICTIONARY}" "${dk2nuTree_HDRS}"
    LINKDEF "${dk2nuTree_LINKDEF}" OPTIONS ${ROOT_DICT_OPT_TREE}
    MODULE dk2nuTree )
endif()

# add dictionary to list of source files
SET(dk2nuTree_SRCS ${dk2nuTree_SRCS} ${dk2nuTree_DICTIONARY})

add_library(dk2nuTree SHARED ${dk2nuTree_SRCS})
target_link_libraries(dk2nuTree ${ROOT_LIBRARIES}
  ${ROOT_Physics_LIBRARY} ${ROOT_Matrix_LIBRARY}
  ${ADD_ROOT_Imt_LIBRARY} ${TBB_LIBRARY} )

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

#file(GLOB dk2nuGenie_HDRS ${PROJECT_SOURCE_DIR}/genie/*.h)
#list(REMOVE_ITEM dk2nuGenie_HDRS ${PROJECT_SOURCE_DIR}/genie/LinkDef.h)
set(dk2nuGenie_HDRS dk2nu/genie/GDk2NuFlux.h)

#MESSAGE("dk2nuGenie_SRCS=${dk2nuGenie_SRCS}")
#MESSAGE(" ")
#MESSAGE("dk2nuGenie_HDRS=${dk2nuGenie_HDRS}")
#MESSAGE(" ")

set(dk2nuGenie_LINKDEF ${PROJECT_SOURCE_DIR}/genie/LinkDef.h)
# do not include ".cxx" here!!!!
set(dk2nuGenie_DICTIONARY dk2nuGenieDict)

#MESSAGE("--DK2NU- about to ROOT_GENERATE_DICTIONARY ${dk2nuGenie_DICTIONARY}")
# option -p = use compiler's preprocessor instead of CINT's
# needed to handle GVersion.h  GRELCODE(2,9,0) construct
if( ${ROOT_MAJOR_VERSION} LESS 6 )
  set(ROOT_DICT_OPT_GENIE -p)
  ROOT_GENERATE_DICTIONARY("${dk2nuGenie_DICTIONARY}" "${dk2nuGenie_HDRS}"
    LINKDEF "${dk2nuGenie_LINKDEF}" OPTIONS ${ROOT_DICT_OPT_GENIE} )
else()
  set(ROOT_INCLUDE_PATH $ENV{ROOT_INCLUDE_PATH}:${PROJECT_SOURCE_DIR}:${PROJECT_SOURCE_DIR}/genie)
  # again define these for "install" phase
  set(GENIEDICTFILES ${CMAKE_BINARY_DIR}/libdk2nuGenie_rdict.pcm
    ${CMAKE_BINARY_DIR}/libdk2nuGenie.rootmap)
  set(ROOT_DICT_OPT_GENIE -I.. -I../genie -I../include -I../include/dk2nu
    -Idk2nu -I../include/genie -I../include/dk2nu/genie -Idk2nu/genie
    -p ) ### -inlineInputHeader )
  # again make use of new MODULE arg
  ROOT_GENERATE_DICTIONARY("${dk2nuGenie_DICTIONARY}" "${dk2nuGenie_HDRS}"
    LINKDEF "${dk2nuGenie_LINKDEF}" OPTIONS ${ROOT_DICT_OPT_GENIE}
    MODULE dk2nuGenie )
endif()

# add dictionary to list of source files
SET(dk2nuGenie_SRCS ${dk2nuGenie_SRCS} ${dk2nuGenie_DICTIONARY})

# don't use -l for dk2nuTree if we want dk2nuGenie to depend on dk2nuTree
# before dk2nuTree is built
add_library(dk2nuGenie SHARED ${dk2nuGenie_SRCS})
target_link_libraries(dk2nuGenie ${ROOT_LIBRARIES}
  ${ROOT_Physics_LIBRARY} ${ROOT_Matrix_LIBRARY}
  dk2nuTree)

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
install(FILES ${TREEDICTFILES} DESTINATION lib)
if(WITH_GENIE)
  install(TARGETS dk2nuGenie DESTINATION lib)
  install(FILES ${GENIEDICTFILES} DESTINATION lib)
endif()
#--------------------
# Install the headers
#   given full dk2nu/tree/*.h names, look up one directory and then down
string(REGEX REPLACE "([^;]+)" "${PROJECT_SOURCE_DIR}/../\\1" dk2nuTree_HDRS2 "${dk2nuTree_HDRS}")
# MESSAGE("dk2nuTree_HDRS2=${dk2nuTree_HDRS2}")
install(FILES   ${dk2nuTree_HDRS2}  DESTINATION include/dk2nu/tree)
if(WITH_GENIE)
  string(REGEX REPLACE "([^;]+)" "${PROJECT_SOURCE_DIR}/../\\1" dk2nuGenie_HDRS2 "${dk2nuGenie_HDRS}")
  # MESSAGE("dk2nuGenie_HDRS2=${dk2nuGenie_HDRS2}")
  install(FILES ${dk2nuGenie_HDRS2} DESTINATION include/dk2nu/genie)
endif()
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
