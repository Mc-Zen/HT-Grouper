#pragma once
#include "graph.h"

namespace Q {

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
}