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

	/// @brief Check if given set of Paulis is measurable with a hardware tailored circuit using one
	///        of the graphs specified. This 
	/// 
	///        NOTE: This version removes all graphs from the given list that cannot be used to diagonalize the Paulis. 
	/// @param collection    Set of pauli operators to diagonalize
	/// @param connectivity  Set of graphs to test for measurability. 
	/// @return success
	template<int n>
	bool is_ht_measurable2(const std::vector<Pauli>& collection, std::vector<Graph<n>>& graphs, HTCircuitFinder& finder) {
		finder.setOperators(collection);
		//println("{}\n{}", collection, connectivity.getAdjacencyMatrix());
		std::vector<int> doErase;
		bool success = false;
		for (const auto& graph : graphs) {
			//println("{}", graph.getAdjacencyMatrix());
			//if (auto result = findHTCircuit(graph, collection)) {
			if (auto result = finder.findHTCircuit(graph)) {
				//println("jo");
				doErase.push_back(false);
				//return true;
				success = true;
				if (collection.size() == 1 || graphs.size() == 1) return true;
			}
			else {
				doErase.push_back(true);
			}
		}
		if (success) {
			std::erase_if(graphs, [&](const auto& graph) { return !finder.findHTCircuit(graph).has_value(); });
		}
		//std::cout << "no\n";
		//println("no, {}", collection);
		return success;
	}


	/// @brief Group Paulis of given hamiltonian into commuting subset that are diagonalizable
	///        with the given hardware connectivity
	/// @param hamiltonian   Hamiltonian specification
	/// @param graphs        Allowed graphs
	/// @return Sets of commuting operators
	template<int n>
	auto applyPauliGrouper(Hamiltonian& hamiltonian, const std::vector<Graph<n>>& graphs) {
		println("Num graphs: {}", graphs.size());

		HTCircuitFinder finder{ hamiltonian.numQubits };

		// Sort by magnitude in descending order 
		std::ranges::sort(hamiltonian.operators, [](const auto& a, const auto& b) {return std::abs(a.second) > std::abs(b.second); });
		//std::ranges::sort(hamiltonian.operators, [](const auto& a, const auto& b) {return std::abs(a.first.pauliWeight()) < std::abs(b.first.pauliWeight()); });


		using Collection = std::pair<std::vector<Pauli>, std::vector<Graph<n>>>;
		std::vector<Collection> collections;
		//std::vector<std::vector<Pauli>> collections;

		int i{};
		for (const auto& [pauli, coefficient] : hamiltonian.operators) {
			bool found{};
			for (Collection& collection : collections) {
				if (!commutesWithAll(collection.first, pauli)) continue;

				collection.first.push_back(pauli);
				if (is_ht_measurable2(collection.first, collection.second, finder)) {
					found = true;
					break;
				}
				collection.first.pop_back();
			}
			if (!found) {
				collections.emplace(collections.end(), std::vector{ pauli }, graphs);
				println("New group ({})", collections.size());
			}
			//std::cout << ++i<<" of " << hamiltonian.operators.size() << "\n";
			println("{} of {}", ++i, hamiltonian.operators.size());
		}
		return collections;
	}
}