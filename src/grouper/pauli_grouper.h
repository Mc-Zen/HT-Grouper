#pragma once

#include <algorithm>
#include "graph.h"
#include "hamiltonian.h"
#include "find_ht_circuit.h"


namespace Q {


	/// @brief Check if given pauli commutes with every other Pauli in the collection. 
	bool commutesWithAll(const std::vector<Pauli>& collection, const Pauli& pauli) {
		for (const auto& p : collection) {
			if (commutator(p, pauli) == 1) return false;
		}
		return true;
	}


	/// @brief Check if given set of Paulis is measurable with a hardware tailored circuit using one
	///        of the graphs specified. 
	/// @param collection    Set of pauli operators to diagonalize
	/// @param connectivity  Set of graphs to test for measurability. 
	/// @return success
	template<int n>
	bool is_ht_measurable(const std::vector<Pauli>& collection, const std::vector<Graph<n>>& graphs) {
		//println("{}\n{}", collection, connectivity.getAdjacencyMatrix());
		for (const auto& graph : graphs) {
			//println("{}", graph.getAdjacencyMatrix());
			if (auto result = findHTCircuit(graph, collection)) {
				println("jo");
				return true;
			}
		}
		println("no, {}", collection);
		return false;
	}


	/// @brief Group Paulis of given hamiltonian into commuting subset that are diagonalizable
	///        with the given hardware connectivity
	/// @param hamiltonian   Hamiltonian specification
	/// @param graphs        Allowed graphs
	/// @return Sets of commuting operators
	template<int n>
	auto applyPauliGrouper(Hamiltonian& hamiltonian, const std::vector<Graph<n>>& graphs) {
		println("Num graphs: {}", graphs.size());
		// Sort by magnitude in descending order 
		std::ranges::sort(hamiltonian.operators, [](const auto& a, const auto& b) {return std::abs(a.second) > std::abs(b.second); });
		//std::ranges::sort(hamiltonian.operators, [](const auto& a, const auto& b) {return std::abs(a.first.pauliWeight()) < std::abs(b.first.pauliWeight()); });

		std::vector<std::vector<Pauli>> collections;

		int i{};
		for (const auto& [pauli, coefficient] : hamiltonian.operators) {
			bool found{};
			for (auto& collection : collections) {
				if (!commutesWithAll(collection, pauli)) continue;

				collection.push_back(pauli);
				if (is_ht_measurable(collection, graphs)) {
					found = true;
					break;
				}
				collection.pop_back();
			}
			if (!found) {
				collections.insert(collections.end(), std::vector{ pauli });
			}
			println("{} of {}", ++i, hamiltonian.operators.size());
		}
		return collections;
	}
}