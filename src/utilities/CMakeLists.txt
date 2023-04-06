
set(target utilities)
add_library(${target}
	string_utility.h
	string_utility.cpp
)
target_include_directories(${target} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


add_unit_test(${target}_unit_tests
	SOURCES 
		tests/string_utility_tests.cpp
	DEPENDENCIES
		${target}
)
