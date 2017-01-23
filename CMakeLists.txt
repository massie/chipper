cmake_minimum_required(VERSION 2.6)
project(chipper)
set (CHIPPER_VERSION_MAJOR 0)
set (CHIPPER_VERSION_MINOR 3)
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

set(LIBLINEAR_ROOT ${CMAKE_SOURCE_DIR}/liblinear-2.1)
include_directories(${LIBLINEAR_ROOT})
set(LIBLINEAR_SRC
    ${LIBLINEAR_ROOT}/linear.cpp
	${LIBLINEAR_ROOT}/linear.h
	${LIBLINEAR_ROOT}/tron.cpp
	${LIBLINEAR_ROOT}/tron.h
	${LIBLINEAR_ROOT}/blas/blas.h
	${LIBLINEAR_ROOT}/blas/daxpy.c
	${LIBLINEAR_ROOT}/blas/ddot.c
	${LIBLINEAR_ROOT}/blas/dnrm2.c
	${LIBLINEAR_ROOT}/blas/dscal.c)
source_group(liblinear FILES ${LIBLINEAR_SRC})

add_subdirectory(src)
include_directories(src)
add_subdirectory(tests)

install(DIRECTORY models/ DESTINATION share/chipper
        PATTERN "models/*"
        PERMISSIONS OWNER_EXECUTE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)

install(FILES tests/test.fa DESTINATION share/chipper
        PERMISSIONS OWNER_EXECUTE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
