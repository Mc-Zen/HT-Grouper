cmake_minimum_required(VERSION 3.19)

project(HT-Grouper)

option(UNIT_TESTING "Enable unit tests for this project" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
find_package(GUROBI REQUIRED)

# Setup gurobi
set(target gurobi_c++)
add_library(${target} INTERFACE)
target_link_libraries(${target} INTERFACE optimized ${GUROBI_CXX_LIBRARY} debug ${GUROBI_CXX_DEBUG_LIBRARY})
target_link_libraries(${target} INTERFACE ${GUROBI_LIBRARY})
target_include_directories(${target} INTERFACE ${GUROBI_INCLUDE_DIRS})

include(cmake/add_unit_test.cmake)

configure_file(cmake/data_path.h.in include/data_path.h)
set(target data)
add_library(${target} INTERFACE)
target_include_directories(${target} INTERFACE "${CMAKE_CURRENT_BINARY_DIR}/include")


if (UNIT_TESTING)
	find_package(Catch2 REQUIRED)
	enable_testing()
	include(Catch)
endif()

add_subdirectory(src)
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT grouper)
