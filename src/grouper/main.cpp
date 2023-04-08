#include "find_ht_circuit.h"
#include "formatting.h"
#include "read_hamiltonians.h"
#include "pauli_grouper.h"

using namespace Q;
using std::cout;



int main() {
	auto filename = R"(C:\Users\alpha\Downloads\hamiltonians.py)"
	auto hamiltonians = readHamiltonians(filename);

	auto& ham = hamiltonians[0];
	auto connectivity = Graph<>::linearGraph(ham.numQubits);

	int maxEdgeCount = 4;
	// Generate all subgraphs of given graph with a maximum of [maxEdgeCount]edges
	auto subgraphs = generateSubgraphs(connectivity, maxEdgeCount); 

	auto collections = applyPauliGrouper(ham, subgraphs);

	println("Found grouping into {} subsets", collections.size());
	for (const auto& collection : collections) {
		println("{}", collection);
	}

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
//	auto connectivityGraph = Graph<>::linearGraph(2);
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