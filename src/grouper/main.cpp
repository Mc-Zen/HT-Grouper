
#include "read_hamiltonians.h"
#include "pauli_grouper.h"
#include "json_formatting.h"
#include "estimated_shot_reduction.h"
#include "data_path.h"
#include "read_config.h"
#include <random>
#include <chrono>

using namespace Q;


std::string toAbsolutePath(const std::string& filename) {
	if (filename.starts_with("C:"))
		return filename;
	return DATA_PATH + filename;
}


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


		using clock = std::chrono::high_resolution_clock;
		const auto t0 = clock::now();

		// Read a hamiltonian consisting of Paulis together with weightings
		// and find a grouping into simultaneously measurable sets respecting
		// a given hardware connectivity. 

		auto filename = toAbsolutePath(config.filename);
		auto outfilename = toAbsolutePath(config.outfilename);
		auto connectivityFile = toAbsolutePath(config.connectivity);

		const auto hamiltonian = readHamiltonianFromJson(filename);
		const auto numQubits = hamiltonian.numQubits;

		Connectivity connectivitySpec = readConnectivity(connectivityFile);
		const auto connectivity = connectivitySpec.getGraph(numQubits);
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

		const auto t1 = clock::now();
		const auto timeInSeconds = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();

		println("Found grouping into {} subsets, run time: {}s", htGrouping.size(), timeInSeconds);


		std::ofstream file{ outfilename };
		auto fileout = std::ostream_iterator<char>(file);

		JsonFormatting::printPauliCollections(fileout, htGrouping, JsonFormatting::MetaInfo{ timeInSeconds, selectedGraphs.size(), connectivity });
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






