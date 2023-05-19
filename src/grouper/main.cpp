
#include "read_hamiltonians.h"
#include "pauli_grouper.h"
#include "python_formatting.h"
#include "json_formatting.h"
#include "estimated_shot_reduction.h"
#include "data_path.h"
#include "read_config.h"
#include <random>


using namespace Q;
using std::cout;



int main() {
	try {

		Configuration config = readConfig(DATA_PATH "config.txt");
		println(R"(Configuration:
  filename = {}
  outfilename = {}
  connectivity = {}
  numThreads = {}
  maxEdgeCount = {}
  numGraphs = {}
  sortGraphsByEdgeCount = {}
)", config.filename, config.outfilename, config.connectivity, config.numThreads, config.maxEdgeCount, config.numGraphs, config.sortGraphsByEdgeCount);



		// Read a hamiltonian consisting of Paulis together with weightings
		// and find a grouping into simultaneously measurable sets respecting
		// a given hardware connectivity. 

		auto filename = config.filename;
		if (!filename.starts_with("C:")) {
			filename = DATA_PATH + filename;
		}
		auto outfilename = config.outfilename;
		if (!outfilename.starts_with("C:")) {
			outfilename = DATA_PATH + outfilename;
		}
		auto connectivityFile = config.connectivity;
		if (!connectivityFile.starts_with("C:")) {
			connectivityFile = DATA_PATH + connectivityFile;
		}

		auto hamiltonian = readHamiltonianFromJson(filename);
		const auto numQubits = hamiltonian.numQubits;

		Connectivity connectivitySpec = readConnectivity(connectivityFile);
		auto connectivity = connectivitySpec.getGraph(numQubits);
		println("Adjacency matrix:\n{}", connectivity.getAdjacencyMatrix());


		// Generate all subgraphs of given graph with a maximum of [maxEdgeCount]edges
		auto subgraphs = generateSubgraphs(connectivity, 0, config.maxEdgeCount);

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
		std::sample(subgraphs.begin(), subgraphs.end(), std::back_inserter(selectedGraphs), config.numGraphs, std::mt19937{ std::random_device{}() });

		if (config.sortGraphsByEdgeCount) {
			std::ranges::sort(selectedGraphs, std::less{}, &Graph<>::edgeCount);
		}

		println("Running pauli grouper with {} Paulis and {} Graphs on {} qubits", hamiltonian.operators.size(), selectedGraphs.size(), numQubits);
		auto htGrouping = applyPauliGrouper2Multithread2(hamiltonian, selectedGraphs, config.numThreads);
		auto tpbGrouping = applyPauliGrouper2Multithread2(hamiltonian, { Graph<>(numQubits) }, config.numThreads, false);

		//htGrouping.erase(htGrouping.begin(), htGrouping.begin() + 2);
		//tpbGrouping.erase(tpbGrouping.begin(), tpbGrouping.begin() + 2);

		auto R_hat_HT = estimated_shot_reduction(hamiltonian, htGrouping);
		auto R_hat_tpb = estimated_shot_reduction(hamiltonian, tpbGrouping);

		println("Found grouping into {} subsets", htGrouping.size());


		std::ofstream file{ outfilename };
		auto out = std::ostream_iterator<char>(std::cout);
		auto fileout = std::ostream_iterator<char>(file);

		//JsonFormatting::printPauliCollections(out, htGrouping);
		JsonFormatting::printPauliCollections(fileout, htGrouping);
		println("Estimated shot reduction\n R_hat_HT = {}\n R_hat_TPB = {}\n R_hat_HT/R_hat_TPB = {}", R_hat_HT, R_hat_tpb, R_hat_HT / R_hat_tpb);
	}
	catch (ConfigReadError& e) {
		println("ConfigReadError: {}", e.what());
	}
	catch (ConnectivityError& e) {
		println("ConnectivityError: {}", e.what());
	}
	catch (std::exception& e) {
		println("{}", e.what());
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