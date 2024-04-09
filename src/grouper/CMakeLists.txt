if (MSVC)
else()
    add_compile_options(-fpermissive)
endif()

set(target grouper)
add_executable(${target} 
	main.cpp
	pauli_grouper.cpp
	read_hamiltonians.h
	pauli_grouper.h
	hamiltonian.h
	python_formatting.h
	json_formatting.h
	estimated_shot_reduction.h
	read_config.h
)
target_link_libraries(${target} PUBLIC q-library gurobi_c++ data)
