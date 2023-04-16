#pragma once

#include <algorithm>
#include "graph.h"
#include "hamiltonian.h"


namespace Q {


	using Collection = std::pair<std::vector<Pauli>, std::vector<Graph<>>>;

	struct CollectionWithGraph {
		std::vector<Pauli> paulis;
		Graph<> graph;
		auto size() const { return paulis.size(); }
	};

	class HTCircuitFinder;


	struct GraphRepr {
		explicit GraphRepr(const Graph<>& graph) : graph(graph), connectedComponents(graph.connectedComponents(true)) {
			for (const auto& component : connectedComponents) {
				uint64_t supportVector{};
				for (auto vertex : component) {
					supportVector |= (1ULL << vertex);
				}
				connectedComponentSupportVectors.push_back(supportVector);
			}
		}

		Graph<> graph;
		std::vector<std::vector<int>> connectedComponents;
		// Support vector for each connected component (a bitstring with 1 
		// for each vertex in the connected component and zeros elsewhere). 
		std::vector<uint64_t> connectedComponentSupportVectors;
	};


	/// @brief Check if given pauli commutes with every other Pauli in the collection. 
	bool commutesWithAll(const std::vector<Pauli>& collection, const Pauli& pauli);

	/// @brief Check if given pauli commutes qubitwise with every other Pauli in the collection. 
	bool qubitwiseCommutesWithAll(const std::vector<Pauli>& collection, const Pauli& pauli);

	/// @brief Check if given pauli commutes locally with every other Pauli in the collection on given support. 
	bool locallyCommutesWithAll(const std::vector<Pauli>& collection, const Pauli& pauli, uint64_t support);

	bool is_ht_measurable(const std::vector<Pauli>& collection, const Graph<>& graph);
	bool is_ht_measurable(const std::vector<Pauli>& collection, const Graph<>& graph) {

	}

	bool is_ht_measurable(const std::vector<Pauli>& collection, const Graph<>& graph, HTCircuitFinder& finder);


	/// @brief Check if given set of Paulis is measurable with a hardware tailored circuit using one
	///        of the graphs specified. 
	/// @param collection    Set of pauli operators to diagonalize
	/// @param connectivity  Set of graphs to test for measurability. 
	/// @return success
	bool is_ht_measurable(const std::vector<Pauli>& collection, const std::vector<Graph<>>& graphs);

	/// @brief Check if given set of Paulis is measurable with a hardware tailored circuit using one
	///        of the graphs specified. This 
	/// 
	///        NOTE: This version removes all graphs from the given list that cannot be used to diagonalize the Paulis. 
	/// @param collection    Set of pauli operators to diagonalize
	/// @param connectivity  Set of graphs to test for measurability. 
	/// @return success
	bool is_ht_measurable2(const std::vector<Pauli>& collection, std::vector<Graph<>>& graphs, HTCircuitFinder& finder);




	/// @brief Group Paulis of given hamiltonian into commuting subset that are diagonalizable
	///        with the given hardware connectivity
	/// @param hamiltonian   Hamiltonian specification
	/// @param graphs        Allowed graphs
	/// @return Sets of commuting operators
	std::vector<Collection> applyPauliGrouper(Hamiltonian& hamiltonian, const std::vector<Graph<>>& graphs);



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
	std::vector<CollectionWithGraph>  applyPauliGrouper2(const Hamiltonian& hamiltonian, const std::vector<Graph<>>& graphs, bool verbose = true);




	/// @brief Group Paulis of given hamiltonian into commuting subset that are diagonalizable
	///		   with the given hardware connectivity. Implementation of Algorithm 1 in https://doi.org/10.48550/arXiv.2203.03646
	/// 
	/// Notes:
	///		- Where two graphs provide results which are considered equally "good" the first graph to occur in 
	///       the graphs argument is selected. Therefore, it may be beneficial to sort these f.e. by edge count. 
	/// 
	/// 
	/// @param hamiltonian   Hamiltonian specification
	/// @param graphs        Allowed graphs (does not need to contains the edgeless graph which will be tested anyway). 
	/// @param verbose       If set to true, will print current status to stdout console output
	/// @return Sets of commuting operators
	std::vector<CollectionWithGraph> applyPauliGrouper2Multithread(const Hamiltonian& hamiltonian, const std::vector<Graph<>>& graphs, int numThreads = 1, bool verbose = true);
}