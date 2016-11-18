# Pass the version info to a header file
configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/project.h.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/project.h")

include(FindZLIB required)

set(CHIPPER_SRC
	aa_props.c
	predict_cleavage.c
	output_fastq.c
	output_netchop.c
	util.c)

source_group(chipper FILES ${CHIPPER_SRC})

add_library(chipper-static STATIC ${CHIPPER_SRC} ${LIBLINEAR_SRC})
target_link_libraries(chipper-static ${ZLIB_LIBRARIES})
set_target_properties(chipper-static PROPERTIES OUTPUT_NAME chipper)

add_executable(chipper chipper.c argtable3.c)
target_link_libraries(chipper chipper-static)

install (TARGETS chipper DESTINATION bin)
install (TARGETS chipper-static DESTINATION lib)
install (FILES chipper.h DESTINATION include)
