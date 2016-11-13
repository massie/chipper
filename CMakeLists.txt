cmake_minimum_required(VERSION 2.6)
project(chipper)
set (CHIPPER_VERSION_MAJOR 0)
set (CHIPPER_VERSION_MINOR 1)
set (CHIPPER_VERSION_PATCH 0)
set (CHIPPER_VERSION "${CHIPPER_VERSION_MAJOR}.${CHIPPER_VERSION_MINOR}.${CHIPPER_VERSION_PATCH}")
enable_testing()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})

include(CPackConfig.txt)
include(CheckCXXCompilerFlag)

#set(CMAKE_VERBOSE_MAKEFILE ON)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

check_cxx_compiler_flag(-Wall temp)
if(temp)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
endif()
check_cxx_compiler_flag(-pedantic temp)
if(temp)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pedantic")
endif()

include_directories(${clipper_SOURCE_DIR}/src)
add_subdirectory(src)
add_subdirectory(tests)
