set(target fmt)

add_library(${target} 
	core.h
	format.h
	format-inl.h
	color.h
	ranges.h
	format.cc
)

target_include_directories(${target} PUBLIC "${CMAKE_CURRENT_LIST_DIR}/..")