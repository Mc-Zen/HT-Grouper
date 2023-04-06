
set(target grouper)
add_executable(${target} 
	main.cpp
	read_hamiltonians.h
	pauli_grouper.h
	hamiltonian.h
)
target_link_libraries(${target} PUBLIC q-library gurobi_c++)
