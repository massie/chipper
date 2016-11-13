# Pass the version info to a header file
configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/version.h.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/version.h")

include(FindZLIB required)

set(CHIPPER_SRC
	aa_props.c
	predict_cleavage.c
	version.h)

source_group(chipper FILES ${CHIPPER_SRC})

add_library(chipper-static STATIC ${CHIPPER_SRC} ${LIBLINEAR_SRC})
target_link_libraries(chipper-static ${ZLIB_LIBRARIES})
set_target_properties(chipper-static PROPERTIES OUTPUT_NAME chipper)

add_executable(chipper chipper.c)
target_link_libraries(chipper chipper-static)

install (TARGETS chipper DESTINATION bin)
