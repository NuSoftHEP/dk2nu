
MESSAGE("--DK2NU-- ONLY_GENIE was false")
set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES}
                        ${PROJECT_SOURCE_DIR}/tree
                        ${PROJECT_SOURCE_DIR}/dk2nu/tree)

include_directories( ${INCLUDE_DIRECTORIES} )
link_directories( ${LINK_DIRECTORIES} )

MESSAGE("--DK2NU-- PROJECT_SOURCE_DIR ${PROJECT_SOURCE_DIR}")
MESSAGE("--DK2NU-- INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES}")


set(LINK_DIRECTORIES ${ROOT_LIBRARY_DIR} ${CMAKE_CURRENT_BINARY_DIR}/lib )
#MESSAGE("CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}")


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
  set(TREEDICTFILES ${CMAKE_BINARY_DIR}/tree/libdk2nuTree_rdict.pcm
                    ${CMAKE_BINARY_DIR}/tree/libdk2nuTree.rootmap )
  set(ROOT_DICT_OPT_TREE -I.. -I../tree
    -I../include -I../include/dk2nu
    -Idk2nu -I../include/tree -I../include/dk2nu/tree -Idk2nu/tree
    -p -inlineInputHeader )
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
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS dk2nuTree DESTINATION lib)
install(FILES ${TREEDICTFILES} DESTINATION lib)

#--------------------
# Install the headers
#   given full dk2nu/tree/*.h names, look up one directory and then down
string(REGEX REPLACE "([^;]+)" "${PROJECT_SOURCE_DIR}/../\\1" dk2nuTree_HDRS2 "${dk2nuTree_HDRS}")
# MESSAGE("dk2nuTree_HDRS2=${dk2nuTree_HDRS2}")
install(FILES   ${dk2nuTree_HDRS2}  DESTINATION include/dk2nu/tree)

