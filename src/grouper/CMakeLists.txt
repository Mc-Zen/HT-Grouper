
set(target grouper)
add_compile_options(-fpermissive)
add_executable(${target} 
	main.cpp
	pauli_grouper.cpp
	read_hamiltonians.h
	random_subgraphs.h
	pauli_grouper.h
	hamiltonian.h
	python_formatting.h
	json_formatting.h
	estimated_shot_reduction.h
	read_config.h
	cli.h
)
target_link_libraries(${target} PUBLIC q-library gurobi_c++ data fmt)
