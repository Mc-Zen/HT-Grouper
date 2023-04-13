#pragma once

#include <algorithm>
#include "graph.h"
#include "hamiltonian.h"
#include "find_ht_circuit.h"
#include <ranges>


namespace Q {


	/// @brief Check if given pauli commutes with every other Pauli in the collection. 
	bool commutesWithAll(const std::vector<Pauli>& collection, const Pauli& pauli) {
		for (const auto& p : collection) {
			if (commutator(p, pauli) == 1) return false;
		}
		return true;
	}

	/// @brief Check if given pauli commutes qubitwise with every other Pauli in the collection. 
	bool qubitwiseCommutesWithAll(const std::vector<Pauli>& collection, const Pauli& pauli) {
		for (const auto& p : collection) {
			if (!commutesQubitWise(p, pauli)) return false;
		}
		return true;
	}

	bool is_ht_measurable(const std::vector<Pauli>& collection, const Graph<>& graph) {
		return findHTCircuit(graph, collection).has_value();
	}

	bool is_ht_measurable(const std::vector<Pauli>& collection, const Graph<>& graph, HTCircuitFinder& finder) {
		finder.setOperators(collection);
		return finder.findHTCircuit(graph).has_value();
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

		bool success = false;
		for (const auto& graph : graphs) {
			if (auto result = finder.findHTCircuit(graph)) {
				success = true;
				if (collection.size() == 1 || graphs.size() == 1) return true;
			}
		}
		if (success) {
			std::erase_if(graphs, [&](const auto& graph) { return !finder.findHTCircuit(graph).has_value(); });
		}
		return success;
	}


	using Collection = std::pair<std::vector<Pauli>, std::vector<Graph<>>>;



	/// @brief Group Paulis of given hamiltonian into commuting subset that are diagonalizable
	///        with the given hardware connectivity
	/// @param hamiltonian   Hamiltonian specification
	/// @param graphs        Allowed graphs
	/// @return Sets of commuting operators
	auto applyPauliGrouper(Hamiltonian& hamiltonian, const std::vector<Graph<>>& graphs) {
		HTCircuitFinder finder{ hamiltonian.numQubits };

		// Sort by magnitude in descending order 
		std::ranges::sort(hamiltonian.operators, [](const auto& a, const auto& b) {return std::abs(a.second) > std::abs(b.second); });

		std::vector<Collection> collections;

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
			if (!found) { // create new group
				collections.emplace(collections.end(), std::vector{ pauli }, graphs);
				println("New group ({})", collections.size());
			}
			println("{} of {}", ++i, hamiltonian.operators.size());
		}
		return collections;
	}

	//std::vector<Pauli>& best(std::vector)

	/// @brief Group Paulis of given hamiltonian into commuting subset that are diagonalizable
	///        with the given hardware connectivity. Implementation of Algorithm 1 in https://doi.org/10.48550/arXiv.2203.03646
	/// Notes:
	///		- Where two graphs provide results which are considered equally "good" the first graph to occur in 
	///       the graphs argument is selected. Therefore, it may be beneficial to sort these f.e. by edge count. 
	/// 
	/// @param hamiltonian   Hamiltonian specification
	/// @param graphs        Allowed graphs (does not need to contains the edgeless graph which will be tested anyway). 
	/// @param verbose       If set to true, will print current status to stdout console output
	/// @return Sets of commuting operators
	auto applyPauliGrouper2(const Hamiltonian& hamiltonian, const std::vector<Graph<>>& graphs, bool verbose = true) {
		HTCircuitFinder finder{ hamiltonian.numQubits };

		auto paulis = hamiltonian.operators;
		// Sort by magnitude in descending order 
		std::ranges::sort(paulis, [](const auto& a, const auto& b) {return std::abs(a.second) > std::abs(b.second); });

		struct Collection {
			std::vector<Pauli> paulis;
			Graph<> graph;
		};
		std::vector<Collection> collections;

		while (!paulis.empty()) {
			const auto& mainPauli = paulis.front().first;
			std::vector<Collection> tempCollections;

			Collection tpbCollection{ { mainPauli }, Graph<>{ hamiltonian.numQubits } };


			for (const auto& [pauli, _] : paulis | std::ranges::views::drop(1)) {
				if (qubitwiseCommutesWithAll(tpbCollection.paulis, pauli)) {
					tpbCollection.paulis.push_back(pauli);
				}
			}
			tempCollections.push_back(tpbCollection);

			int j{};
			for (const auto& graph : graphs) {
				if (verbose) print("\33[2K\rGraph {:>4} of {:>4}", j++, graphs.size());

				Collection collection{ {mainPauli }, graph };
				if (!is_ht_measurable(collection.paulis, graph, finder)) continue;

				for (const auto& [pauli, _] : paulis | std::ranges::views::drop(1)) {
					if (!commutesWithAll(collection.paulis, pauli)) continue;
					collection.paulis.push_back(pauli);
					if (!is_ht_measurable(collection.paulis, graph, finder)) {
						collection.paulis.pop_back();
					}
				}
				tempCollections.push_back(collection);
			}

			const auto& bestCollection = *std::ranges::max_element(tempCollections, std::less{}, [](auto& c) {return c.paulis.size(); });
			collections.push_back(bestCollection);
			for (const auto& pauli : bestCollection.paulis) {
				std::erase_if(paulis, [&pauli](auto& val) {return val.first == pauli; });
			}
			if (verbose) println("\33[2K\r{} of {} remaining ({} group{})", paulis.size(), hamiltonian.operators.size(), collections.size(), collections.size() == 1 ? "" : "s");
		}
		return collections;
	}
}