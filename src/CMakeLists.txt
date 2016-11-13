# Pass the version info to a header file
configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/version.h.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/version.h")

include(FindZLIB required)

set(CHIPPER_SRC
	aa_props.c
	predict_cleavage.c
	version.h)

set(LIBLINEAR_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/liblinear-2.1)
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

source_group(chipper FILES ${CHIPPER_SRC})
source_group(liblinear FILES ${LIBLINEAR_SRC})

include_directories(${LIBLINEAR_ROOT})
add_library(chipper-static STATIC ${CHIPPER_SRC} ${LIBLINEAR_SRC})
target_link_libraries(chipper-static ${ZLIB_LIBRARIES})
set_target_properties(chipper-static PROPERTIES OUTPUT_NAME chipper)

add_executable(chipper chipper.c)
target_link_libraries(chipper chipper-static)

install (TARGETS chipper DESTINATION bin)
