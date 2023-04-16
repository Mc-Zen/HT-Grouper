
set(target grouper)
add_executable(${target} 
	main.cpp
	pauli_grouper.cpp
	read_hamiltonians.h
	pauli_grouper.h
	hamiltonian.h
	python_formatting.h
)
target_link_libraries(${target} PUBLIC q-library gurobi_c++)
