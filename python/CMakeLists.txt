cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

find_package(Python3 REQUIRED COMPONENTS Development Interpreter)

# PYTHON PATHS
set(DK2NU_PYTHONPATH "python/${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}")

include(FetchContent)

FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        f5fbe86 # release-3.0.1
)

FetchContent_MakeAvailable(pybind11)

pybind11_add_module(pydk2nu SHARED pydk2nu.cxx)

target_include_directories(pydk2nu
  PRIVATE 
    ${CMAKE_CURRENT_LIST_DIR}/../tree)

target_link_libraries(pydk2nu 
  PRIVATE 
    dk2nuTree)

install(TARGETS pydk2nu
    LIBRARY DESTINATION ${DK2NU_PYTHONPATH})