macro(add_chipper_executable name)
    set(source "${ARGV1}")
    if (NOT source)
        set(source "${name}.c")
    endif (NOT source)
    add_executable(${name} ${source})
    target_link_libraries(${name} chipper-static)
endmacro(add_chipper_executable)

macro(add_chipper_test name)
    add_chipper_executable(${name} ${ARGN})
	set(exec_name ${EXECUTABLE_OUTPUT_PATH}/${name})
    add_test(${name}
        ${CMAKE_COMMAND} -E chdir ${CMAKE_SOURCE_DIR}/tests
        ${exec_name}
    )
endmacro(add_chipper_test)

include_directories(${CMAKE_SOURCE_DIR}/src)

add_chipper_test(test_predict_cleavage)
add_chipper_test(test_aa_props)
