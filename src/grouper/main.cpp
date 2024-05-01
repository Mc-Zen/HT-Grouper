
#include "read_hamiltonians.h"
#include "pauli_grouper.h"
#include "json_formatting.h"
#include "estimated_shot_reduction.h"
#include "data_path.h"
#include "read_config.h"
#include "random_subgraphs.h"
#include <random>
#include <chrono>
#include <filesystem>
#include <iostream>

using namespace Q;


std::string toAbsolutePath(const std::string& filename) {
	if (filename.starts_with("C:") || filename.starts_with("/"))
		return filename;
	return DATA_PATH + filename;
}




int main(int argc, char** argv) {
	try {
		std::string configPath = DATA_PATH "config.txt";
		if (argc == 2) {
			configPath = argv[1];
		}
		Configuration config = readConfig(configPath);
		fmt::println(R"(Configuration:
  filename = {}
  outfilename = {}
  connectivity = {}
  numThreads = {}
  maxEdgeCount = {}
  numGraphs = {}
  sortGraphsByEdgeCount = {}
  extractComputationalBasis = {}
	generateTPBs = {}
)", config.filename, config.outfilename, config.connectivity, config.numThreads, config.maxEdgeCount, config.numGraphs, config.sortGraphsByEdgeCount,
config.extractComputationalBasis, config.generateTPBs);


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
		fmt::println("Adjacency matrix:\n{}", connectivity.getAdjacencyMatrix());


		std::cout << std::endl;

		auto outPath = std::filesystem::path(outfilename);
		std::filesystem::create_directories(outPath.parent_path());

		// Generate all subgraphs of given graph with a maximum of [maxEdgeCount]edges
		//auto subgraphs = generateSubgraphs(connectivity, 0, config.maxEdgeCount);

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
		//	fmt::println("{}", subgraphs[i-128].getAdjacencyMatrix());
		//}



		const auto seed = config.seed == 0 ? std::random_device{}() : config.seed;
		std::mt19937_64 randomGenerator{ seed };

		JsonFormatting::MetaInfo metaInfo{
			.randomSeed = seed,
			.connectivity = connectivity,
			.inputFilename = std::filesystem::path(config.filename).filename().string(),
			.grouperType = config.grouperType
		};

		//decltype(subgraphs) selectedGraphs;
		//std::sample(subgraphs.begin(), subgraphs.end(), std::back_inserter(selectedGraphs), config.numGraphs, randomGenerator);
		auto selectedGraphs = getRandomSubgraphs(connectivity, config.numGraphs, config.maxEdgeCount, randomGenerator);

		if (config.sortGraphsByEdgeCount) {
			std::ranges::sort(selectedGraphs, std::less{}, &Graph<>::edgeCount);
		}

		fmt::println("Running HT Pauli grouper with {} Paulis and {} Graphs on {} qubits", hamiltonian.operators.size(), selectedGraphs.size(), numQubits);
		fmt::println("Random seed: {}\n", seed);

		std::unique_ptr<PauliGrouper> grouper;
		if (config.grouperType == 1) {
			grouper = std::make_unique<PauliGrouper>(hamiltonian, selectedGraphs, config.numThreads, config.extractComputationalBasis, config.verboseLog);
		}
		else if (config.grouperType == 2) {
			grouper = std::make_unique<PauliGrouper2>(hamiltonian, connectivity, config.numThreads, config.extractComputationalBasis, config.verboseLog, seed, config.numGraphs);
			//PauliGrouper grouper(hamiltonian, selectedGraphs, config.numThreads, config.extractComputationalBasis, config.verboseLog);
		}
		else {
			throw std::runtime_error(fmt::format("Really bad hazard warning. "));
		}
		int count{};
		while (*grouper) {
			++count;
			auto collection = grouper->groupOne();
			if (config.intermediateFileFrequency != 0 && count % config.intermediateFileFrequency == 0) {

				auto parentPath = outPath.parent_path().string();
				auto stem = outPath.stem().string();
				auto ext = outPath.extension().string();

				std::ofstream file{ parentPath + "/" + stem + "_savingpoint" + ext };
				auto fileout = std::ostream_iterator<char>(file);

				const auto tTemp2 = clock::now();
				const auto timeInSeconds = std::chrono::duration_cast<std::chrono::seconds>(tTemp2 - t0).count();
				auto tempMetaInfo = metaInfo;
				tempMetaInfo.timeInSeconds = timeInSeconds;
				tempMetaInfo.numGraphs = selectedGraphs.size();
				JsonFormatting::printPauliCollections(fileout, grouper->getCollections(), tempMetaInfo);
			}
		}

		const auto& htGrouping = grouper->getCollections();
		//auto htGrouping = applyPauliGrouper2Multithread2(hamiltonian, selectedGraphs, config.numThreads, config.extractComputationalBasis);


		auto R_hat_HT = estimated_shot_reduction(hamiltonian, htGrouping);

		double R_hat_tpb = 0;
		if (config.generateTPBs) {
			fmt::println("\n\n\n---------------\nRunning TPB grouping");
			auto tpbGrouping = applyPauliGrouper2Multithread2(hamiltonian, { Graph<>(numQubits) }, config.numThreads, true, config.verboseLog);
			R_hat_tpb = estimated_shot_reduction(hamiltonian, tpbGrouping);
			fmt::println("\n\n\n---------------\n");
		}

		//htGrouping.erase(htGrouping.begin(), htGrouping.begin() + 2);
		//tpbGrouping.erase(tpbGrouping.begin(), tpbGrouping.begin() + 2);


		const auto t1 = clock::now();
		const auto timeInSeconds = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();

		fmt::println("Found grouping into {} subsets, run time: {}s", htGrouping.size(), timeInSeconds);


		std::ofstream file{ outPath };
		auto fileout = std::ostream_iterator<char>(file);

		metaInfo.timeInSeconds = timeInSeconds;
		metaInfo.numGraphs = selectedGraphs.size();
		metaInfo.Rhat_HT = R_hat_HT;
		metaInfo.Rhat_TPB = R_hat_tpb;
		JsonFormatting::printPauliCollections(fileout, htGrouping, metaInfo);
		fmt::println("Estimated shot reduction\n R_hat_HT = {}\n R_hat_TPB = {}\n R_hat_HT/R_hat_TPB = {}", R_hat_HT, R_hat_tpb, R_hat_HT / R_hat_tpb);
		std::cout << std::endl;
	}
	catch (ConfigReadError& e) {
		fmt::println("ConfigReadError: {}", e.what());
	}
	catch (ConnectivityError& e) {
		fmt::println("ConnectivityError: {}", e.what());
	}
	catch (std::exception& e) {
		fmt::println("{}", e.what());
	}
	return 0;
}






