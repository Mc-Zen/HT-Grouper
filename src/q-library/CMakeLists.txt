set(target q-library)
add_library(${target}
	basic_operators.h
	binary.h
	binary_pauli.h
	binary_phase.h
	efficient_mub.h
	evolve_pauli.h
	find_ht_circuit.h
	formatting.h
	graph.h
	ht_circuits.h
	lc_classes.h
	matrix.h
	n_choose_2_iterator.h
	pauli.h
	pauli_operator_map.h
	sector_length_distribution.h
	special_math.h
	stabilizer.h
	symbolic.h
	symbolic.cpp
	random_clifford_element.h
	mub_analysis.h
	get_H_base_from_irreducible_polynomial.h
	generate_mub.h
	quantum_circuit.h
)
target_include_directories(${target} PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(${target} PUBLIC utilities)

add_unit_test(${target}_unit_tests
	SOURCES 
		tests/graph_tests.cpp
		tests/sector_length_distribution_tests.cpp
		tests/efficient_binary_math_tests.cpp
		tests/binary_pauli_tests.cpp
		tests/lc_classes_tests.cpp
		tests/matrix_tests.cpp
		tests/pauli_tests.cpp
	DEPENDENCIES
		${target}
)
