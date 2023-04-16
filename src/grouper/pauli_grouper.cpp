#include "pauli_grouper.h"

#include "find_ht_circuit.h"
#include <ranges>
#include <thread>

using namespace Q;


/// @brief Check if given pauli commutes with every other Pauli in the collection. 
bool Q::commutesWithAll(const std::vector<Pauli>& collection, const Pauli& pauli) {
	for (const auto& p : collection) {
		if (commutator(p, pauli) == 1) return false;
	}
	return true;
}

/// @brief Check if given pauli commutes qubitwise with every other Pauli in the collection. 
bool Q::qubitwiseCommutesWithAll(const std::vector<Pauli>& collection, const Pauli& pauli) {
	for (const auto& p : collection) {
		if (!commutesQubitWise(p, pauli)) return false;
	}
	return true;
}

bool Q::locallyCommutesWithAll(const std::vector<Pauli>& collection, const Pauli& pauli, uint64_t support)
{
	for (const auto& p : collection) {
		if (!commutesLocally(p, pauli, support)) return false;
	}
	return true;
}

bool Q::is_ht_measurable(const std::vector<Pauli>& collection, const Graph<>& graph) {
	return findHTCircuit(graph, collection).has_value();
}

bool Q::is_ht_measurable(const std::vector<Pauli>& collection, const Graph<>& graph, HTCircuitFinder& finder) {
	finder.setOperators(collection);
	return finder.findHTCircuit(graph).has_value();
}

/// @brief Check if given set of Paulis is measurable with a hardware tailored circuit using one
///        of the graphs specified. 
/// @param collection    Set of pauli operators to diagonalize
/// @param connectivity  Set of graphs to test for measurability. 
/// @return success
bool Q::is_ht_measurable(const std::vector<Pauli>& collection, const std::vector<Graph<>>& graphs) {
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
bool Q::is_ht_measurable2(const std::vector<Pauli>& collection, std::vector<Graph<>>& graphs, HTCircuitFinder& finder) {
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

/// @brief Group Paulis of given hamiltonian into commuting subset that are diagonalizable
///        with the given hardware connectivity
/// @param hamiltonian   Hamiltonian specification
/// @param graphs        Allowed graphs
/// @return Sets of commuting operators
std::vector<Collection> Q::applyPauliGrouper(Hamiltonian& hamiltonian, const std::vector<Graph<>>& graphs) {
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
std::vector<CollectionWithGraph> Q::applyPauliGrouper2(const Hamiltonian& hamiltonian, const std::vector<Graph<>>& graphs, bool verbose) {
	HTCircuitFinder finder{ hamiltonian.numQubits };

	auto paulis = hamiltonian.operators;
	// Sort by magnitude in descending order 
	std::ranges::sort(paulis, [](const auto& a, const auto& b) {return std::abs(a.second) > std::abs(b.second); });

	std::vector<CollectionWithGraph> collections;

	while (!paulis.empty()) {
		const auto& mainPauli = paulis.front().first;
		std::vector<CollectionWithGraph> tempCollections;

		CollectionWithGraph tpbCollection{ { mainPauli }, Graph<>{ hamiltonian.numQubits } };


		for (const auto& [pauli, _] : paulis | std::ranges::views::drop(1)) {
			if (qubitwiseCommutesWithAll(tpbCollection.paulis, pauli)) {
				tpbCollection.paulis.push_back(pauli);
			}
		}
		tempCollections.push_back(tpbCollection);

		int j{};
		for (const auto& graph : graphs) {
			if (verbose) print("\33[2K\rGraph {:>4} of {:>4}", j++, graphs.size());

			CollectionWithGraph collection{ { mainPauli }, graph };
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
std::vector<CollectionWithGraph> Q::applyPauliGrouper2Multithread(const Hamiltonian& hamiltonian, const std::vector<Graph<>>& graphs, int numThreads, bool verbose) {
	const auto numGraphsPerThread = static_cast<size_t>(std::ceil(static_cast<float>(graphs.size()) / numThreads));
	std::vector<HTCircuitFinder> finders;
	for (int i = 0; i < numThreads; ++i) finders.emplace_back(hamiltonian.numQubits);

	auto paulis = hamiltonian.operators;
	// Sort by magnitude in descending order 
	std::ranges::sort(paulis, [](const auto& a, const auto& b) {return std::abs(a.second) > std::abs(b.second); });

	std::vector<CollectionWithGraph> collections;

	while (!paulis.empty()) {
		const auto& mainPauli = paulis.front().first;

		CollectionWithGraph tpbCollection{ { mainPauli }, Graph<>{ hamiltonian.numQubits } };

		for (const auto& [pauli, _] : paulis | std::ranges::views::drop(1)) {
			if (qubitwiseCommutesWithAll(tpbCollection.paulis, pauli)) {
				tpbCollection.paulis.push_back(pauli);
			}
		}

		std::atomic_int visitedGraphs{};
		std::atomic_int finishedThreads{};

		auto work = [&](size_t first, size_t last, std::vector<CollectionWithGraph>& partialSolution, HTCircuitFinder& finder) {
			for (auto i = first; i < last; ++i) {
				++visitedGraphs;
				const auto& graph = graphs[i];
				CollectionWithGraph collection{ { mainPauli }, graph };
				if (!is_ht_measurable(collection.paulis, graph, finder)) continue;

				for (const auto& [pauli, _] : paulis | std::ranges::views::drop(1)) {
					if (!commutesWithAll(collection.paulis, pauli)) continue;
					collection.paulis.push_back(pauli);
					if (!is_ht_measurable(collection.paulis, graph, finder)) {
						collection.paulis.pop_back();
					}
				}
				partialSolution.push_back(collection);
			}
			++finishedThreads;
		};

		std::vector<std::vector<CollectionWithGraph>> partialSolutions(numThreads);

		{
			std::vector<std::jthread> workers;
			for (int i = 0; i < numThreads; ++i) {
				const auto firstGraphIndex = numGraphsPerThread * i;
				const auto lastGraphIndex = numGraphsPerThread * (i + 1);
				workers.emplace_back(work, firstGraphIndex, std::min(lastGraphIndex, graphs.size()), std::ref(partialSolutions[i]), std::ref(finders[i]));
			}

			if (verbose) {
				int previousVisitedGraphs = -1;
				while (finishedThreads < numThreads) {
					if (int currentlyVisitedGraphs = visitedGraphs.load(); currentlyVisitedGraphs != previousVisitedGraphs) {
						print("\33[2K\rGraph {:>4} of {:>4}", visitedGraphs.load(), graphs.size());
						previousVisitedGraphs = currentlyVisitedGraphs;
					}
					using namespace std::chrono_literals;
					std::this_thread::sleep_for(10ms);
				}
			}
		}

		const auto* bestCollection = &tpbCollection;
		for (const auto& partialSolution : partialSolutions) {
			for (const auto& collection : partialSolution) {
				if (collection.size() > bestCollection->size()) bestCollection = &collection;
			}
		}
		collections.push_back(*bestCollection);
		for (const auto& pauli : bestCollection->paulis) {
			std::erase_if(paulis, [&pauli](auto& val) { return val.first == pauli; });
		}
		if (verbose) println("\33[2K\r{} of {} remaining ({} group{}): {} -> {}\n",
			paulis.size(), hamiltonian.operators.size(), collections.size(), collections.size() == 1 ? "" : "s",
			collections.back().paulis, collections.back().graph.getEdges());
	}
	return collections;
}
