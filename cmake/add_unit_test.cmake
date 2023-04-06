
# Add Catch2 unit tests by creating a target. 
#
# This function links the created target against Catch2WithMain. 
# The created test target is automatically "discovered" by calling catch_discover_tests(). 
# ____________________________________
#
# Parameters:
#   [target]        target name
#   SOURCES	        input source .cpp files
#   DEPENDENCIES    targets to link against
#   FOLDER          IDE folder to set, defaults to "Unit Tests" if not provided
# ____________________________________
#
# Call example: 
# 
#	 add_unit_test(synth_unit_tests
#		SOURCES 
#			tests/wavetable_tests.cpp
#			tests/wavetable_oscillator_tests.cpp
#		DEPENDENCIES
#			synth
#		FOLDER
#			"My Unit Tests"
#	 )

function(add_unit_test target)
	set(oneValueArgs FOLDER)
	set(multiValueArgs SOURCES DEPENDENCIES)
	cmake_parse_arguments(PARSE_ARGV 0 PARAMS "${options}" "${oneValueArgs}" "${multiValueArgs}")

	
	if(UNIT_TESTING)
		if(NOT PARAMS_FOLDER)
			set(PARAMS_FOLDER "Unit Tests")
		endif()

		add_executable(${target} ${PARAMS_SOURCES})
		target_link_libraries(${target} PRIVATE Catch2::Catch2WithMain)
		target_link_libraries(${target} PRIVATE ${PARAMS_DEPENDENCIES})
		set_target_properties(${target} PROPERTIES FOLDER ${PARAMS_FOLDER})
		catch_discover_tests(${target})
	endif()

endfunction()
