#include "find_ht_circuit.h"
#include "formatting.h"
#include "read_hamiltonians.h"
#include "pauli_grouper.h"
#include <random>


using namespace Q;
using std::cout;



int main() {

	//HTCircuitFinder finder{ 12 , true};
	//auto result = is_ht_measurable(std::vector<Pauli>{Pauli{ "IYXIIIIIXZIX" }}, Graph<>{12}, finder);

	//auto filename = R"(C:\Users\E-Bow\Documents\Code\Cplusplus\HT-Grouper\data\hamiltonians.py)";
	auto filename = R"(C:\Users\alpha\Downloads\hamiltonians.py)";
	auto hamiltonians = readHamiltonians(filename);

	auto& ham = hamiltonians[2];
	auto connectivity = Graph<>::linear(ham.numQubits);

	int maxEdgeCount = 90;
	// Generate all subgraphs of given graph with a maximum of [maxEdgeCount]edges
	auto subgraphs = generateSubgraphs(connectivity, 0, maxEdgeCount);
	const auto numQubits = ham.numQubits;

	//for (int starSize = 3; starSize < 8; ++starSize) {
	//	for (int start = 0; start <= numQubits - starSize; ++start) {
	//		Graph<> star{ numQubits };
	//		for (int i = 0; i < starSize; ++i) {
	//			star.addEdge(start, start + i);
	//		}
	//		subgraphs.push_back(star);
	//	}
	//}

	//std::ranges::rotate(subgraphs, subgraphs.begin() + 128);
	//for (int i = 128; i < subgraphs.size(); ++i) {
	//	println("{}", subgraphs[i-128].getAdjacencyMatrix());
	//}


	decltype(subgraphs) selectedGraphs;
	std::sample(subgraphs.begin(), subgraphs.end(), std::back_inserter(selectedGraphs), 50, std::mt19937{ std::random_device{}() });
	std::ranges::sort(selectedGraphs, std::less{}, &Graph<>::edgeCount);

	for (auto g : selectedGraphs)println("{}", g.edgeCount());
	println("Running pauli grouper with {} Paulis and {} Graphs on {} qubits", ham.operators.size(), selectedGraphs.size(), numQubits);
	auto collections = applyPauliGrouper2(ham, selectedGraphs);

	//println("Found grouping into {} subsets", collections.size());
	//for (const auto& collection : collections) {
	//	println("{}", collection.first);
	//	println("{} ({} graphs)", collection.second[0].getAdjacencyMatrix(), collection.second.size());
	//}
	println("Found grouping into {} subsets", collections.size());
	for (const auto& collection : collections) {
		println("{} -> {}", collection.paulis, collection.graph.getEdges());
	}

	return 0;
}






int main2() {
	auto filename = R"(C:\Users\E-Bow\Downloads\grouping H2O2.txt)";
	auto groups = readPauliGroups(filename);

	std::vector<Pauli> singles;

	int initialGroupCount{};
	for (const auto& group : groups) {
		if (group.size() <= 2) {
			std::ranges::copy(group, std::back_inserter(singles));

			++initialGroupCount;
		}
		//println("{}", group);
	}
	auto numQubits = singles[0].numQubits();

	auto connectivity = Graph<>::linear(numQubits);

	int maxEdgeCount = 7;
	// Generate all subgraphs of given graph with a maximum of [maxEdgeCount]edges
	auto subgraphs = generateSubgraphs(connectivity, 1, maxEdgeCount);
	decltype(subgraphs) selectedGraphs;
	selectedGraphs.push_back(Graph<>(numQubits));
	std::sample(subgraphs.begin(), subgraphs.end(), std::back_inserter(selectedGraphs), 10, std::mt19937{ std::random_device{}() });

	Hamiltonian ham;
	ham.numQubits = numQubits;
	for (const auto& pauli : singles) {
		ham.operators.emplace_back(pauli, 0.);
	}

	println("trying to regroup {} groups with {} Paulis", initialGroupCount, singles.size());
	auto collections = applyPauliGrouper(ham, selectedGraphs);
	return 0;
}




// some testing
//void findReadoutCircuit2Qubits() {
//	constexpr auto set1 = parseMubSet<2, 3>("XY YZ ZX");
//	constexpr auto set2 = parseMubSet<2, 3>("XX YY ZZ");
//	constexpr auto set3 = parseMubSet<2, 3>("XI XZ IZ");
//	constexpr auto set4 = parseMubSet<2, 3>("YI YX IX");
//	constexpr auto set5 = parseMubSet<2, 3>("ZI ZY IY");
//	const std::vector a = { &set1,&set2,&set3,&set4,&set5 };
//
//
//	auto connectivityGraph = Graph<>::linear(2);
//
//	const auto graphs = generateSubgraphs(connectivityGraph);
//	for (const auto& base : a) {
//		for (const auto& graph : graphs) {
//			const auto result = findHTCircuit(graph, *base);
//			if (result) {
//				const auto& layer = result.value();
//				println("Circuit found\n{}", graph.getAdjacencyMatrix());
//				break;
//			}
//		}
//	}
//}