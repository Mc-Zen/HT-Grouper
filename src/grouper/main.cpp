
#include "read_hamiltonians.h"
#include "pauli_grouper.h"
#include "json_formatting.h"
#include "estimated_shot_reduction.h"
#include "data_path.h"
#include "read_config.h"
#include <random>
#include <chrono>
#include <filesystem>

using namespace Q;


std::string toAbsolutePath(const std::string& filename) {
	if (filename.starts_with("C:") || filename.starts_with("/"))
		return filename;
	return DATA_PATH + filename;
}


template<class RNG>
auto getRandomSubgraphs(const Graph<>& graph, int64_t num, int maxEdgeCount, RNG&& rng) {
	const auto edgeCount = graph.edgeCount();
	if (edgeCount <= 63) {
		// Check if num wanted graphs is greater or equal the total number of subgraphs
		// then we just return all subgraphs 
		const uint64_t totalNumSubgraphs = 1ULL << edgeCount;
		if (num >= totalNumSubgraphs) {
			return generateSubgraphs(graph, 0, maxEdgeCount);
		}
		else {
			auto edges = graph.getEdges();
			auto edgeMask = (1ULL << edgeCount) - 1;
			std::vector<Graph<>> subgraphs;
			while (subgraphs.size() < num) {
				uint64_t randomInt = rng() & edgeMask;
				if (const auto ec = std::popcount(randomInt); ec > maxEdgeCount) continue;

				Graph<> subgraph(graph.graphSize);
				for (size_t j = 0; j < edges.size(); ++j) {
					if (randomInt & (1ULL << j)) {
						subgraph.addEdge(edges[j].first, edges[j].second);
					}
				}
				subgraphs.push_back(subgraph);
			}
			return subgraphs;
		}
	}
	else {
		throw std::runtime_error("More than 63 edges are currently not supported");
	}
}


int main(int argc, char**argv) {
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
		fmt::println("Adjacency matrix:\n{}", connectivity.getAdjacencyMatrix());



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
		//decltype(subgraphs) selectedGraphs;
		//std::sample(subgraphs.begin(), subgraphs.end(), std::back_inserter(selectedGraphs), config.numGraphs, randomGenerator);
		auto selectedGraphs = getRandomSubgraphs(connectivity, config.numGraphs, config.maxEdgeCount, randomGenerator);

		if (config.sortGraphsByEdgeCount) {
			std::ranges::sort(selectedGraphs, std::less{}, &Graph<>::edgeCount);
		}

		fmt::println("Running HT Pauli grouper with {} Paulis and {} Graphs on {} qubits", hamiltonian.operators.size(), selectedGraphs.size(), numQubits);
		fmt::println("Random seed: {}\n", seed);

		PauliGrouper grouper(hamiltonian, selectedGraphs, config.numThreads, config.extractComputationalBasis);
		int count{};
		while (grouper) {
			++count;
			auto collection = grouper.groupOne();
			if (config.intermediateFileFrequency != 0 && count % config.intermediateFileFrequency == 0) {

				auto parentPath = outPath.parent_path().string();
				auto stem = outPath.stem().string();
				auto ext = outPath.extension().string();

				std::ofstream file{ parentPath + "/" + stem + "_savingpoint_" + std::to_string(count) + ext };
				auto fileout = std::ostream_iterator<char>(file);

				const auto tTemp2 = clock::now();
				const auto timeInSeconds = std::chrono::duration_cast<std::chrono::seconds>(tTemp2 - t0).count();
				JsonFormatting::printPauliCollections(fileout, grouper.getCollections(), JsonFormatting::MetaInfo{ timeInSeconds, selectedGraphs.size(), seed, connectivity });
			}
		}

		const auto& htGrouping = grouper.getCollections();
		//auto htGrouping = applyPauliGrouper2Multithread2(hamiltonian, selectedGraphs, config.numThreads, config.extractComputationalBasis);


		auto R_hat_HT = estimated_shot_reduction(hamiltonian, htGrouping);

		double R_hat_tpb = 0;
		if (config.generateTPBs) {
			fmt::println("\n\n\n---------------\nRunning TPB grouping", hamiltonian.operators.size(), selectedGraphs.size(), numQubits);
			auto tpbGrouping = applyPauliGrouper2Multithread2(hamiltonian, { Graph<>(numQubits) }, config.numThreads, false);
			R_hat_tpb = estimated_shot_reduction(hamiltonian, tpbGrouping);
		}

		//htGrouping.erase(htGrouping.begin(), htGrouping.begin() + 2);
		//tpbGrouping.erase(tpbGrouping.begin(), tpbGrouping.begin() + 2);


		const auto t1 = clock::now();
		const auto timeInSeconds = std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();

		fmt::println("Found grouping into {} subsets, run time: {}s", htGrouping.size(), timeInSeconds);


		std::ofstream file{ outPath };
		auto fileout = std::ostream_iterator<char>(file);

		JsonFormatting::printPauliCollections(fileout, htGrouping, JsonFormatting::MetaInfo{
			.timeInSeconds = timeInSeconds,
			.numGraphs = selectedGraphs.size(),
			.randomSeed = seed,
			.connectivity = connectivity,
			.Rhat_HT = R_hat_HT,
			.Rhat_TPB = R_hat_tpb,
			});
		fmt::println("Estimated shot reduction\n R_hat_HT = {}\n R_hat_TPB = {}\n R_hat_HT/R_hat_TPB = {}", R_hat_HT, R_hat_tpb, R_hat_HT / R_hat_tpb);
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






