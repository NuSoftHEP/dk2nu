
set(GENIE $ENV{GENIE})
MESSAGE("--DK2NU-- WITH_GENIE or GENIE_ONLY was true")
if (GENIE_ONLY)
  if ( "${DK2NUDATA_DIR}" STREQUAL "" )
    MESSAGE("DK2NUDATA_DIR is NULL, GENIE_ONLY is set ...")
    MESSAGE("this will likely end in tears")
    MESSAGE("===============================================================")
    endif()
  set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES}
                          ${DK2NUDATA_DIR}/include)
else()
  set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES}
                          ${PROJECT_SOURCE_DIR}/tree
                          ${PROJECT_SOURCE_DIR}/dk2nu/tree)
endif()

set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES}
                        ${PROJECT_SOURCE_DIR}/genie
                        ${PROJECT_SOURCE_DIR}/dk2nu/genie
                        $ENV{LIBXML2_INC}
                        $ENV{LOG4CPP_INC}
                        ${GENIE}/src)

set(LINK_DIRECTORIES ${LINK_DIRECTORIES} ${GENIE}/lib )
if (GENIE_ONLY)
  set(LINK_DIRECTORIES ${LINK_DIRECTORIES} ${DK2NUDATA_DIR}/lib )
    MESSAGE("DK2NUDATA_DIR = ${DK2NUDATA_DIR}")
    MESSAGE("adding to LINK_DIRECTORIES => ${LINK_DIRECTORIES}")
endif()

include_directories( ${INCLUDE_DIRECTORIES} )
link_directories( ${LINK_DIRECTORIES} )


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
  set(GENIEDICTFILES ${CMAKE_BINARY_DIR}/genie/libdk2nuGenie_rdict.pcm
                     ${CMAKE_BINARY_DIR}/genie/libdk2nuGenie.rootmap )
  set(ROOT_DICT_OPT_GENIE -I.. -I../genie -I../include -I../include/dk2nu
    -Idk2nu -I../include/genie -I../include/dk2nu/genie -Idk2nu/genie
    -p -inlineInputHeader )
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
  ${ADD_ROOT_Imt_LIBRARY} ${TBB_LIBRARY}
  dk2nuTree)

#set_property(TARGET dk2nuGenie PROPERTY DEPENDS dk2nuTree)
#add_custom_command(OUTPUT dk2nuGenie COMMAND echo hey there DEPENDS dk2nuTree)

#MESSAGE("--DK2NU- dk2nuGenie section done")

install(TARGETS dk2nuGenie DESTINATION lib)
install(FILES ${GENIEDICTFILES} DESTINATION lib)


string(REGEX REPLACE "([^;]+)" "${PROJECT_SOURCE_DIR}/../\\1" dk2nuGenie_HDRS2 "${dk2nuGenie_HDRS}")
# MESSAGE("dk2nuGenie_HDRS2=${dk2nuGenie_HDRS2}")
install(FILES ${dk2nuGenie_HDRS2} DESTINATION include/dk2nu/genie)
