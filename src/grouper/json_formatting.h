#pragma once
#include "formatting.h"
#include "graph.h"

namespace JsonFormatting {


	struct MetaInfo {
		long long timeInSeconds{};
		size_t numGraphs{};
		size_t randomSeed{};
		Q::Graph<> connectivity;
		double Rhat_HT{};
		double Rhat_TPB{};
	};

	void printEdgeList(auto out, const std::vector<std::pair<int, int>>& edges) {
		for (size_t i = 0; i < edges.size(); ++i) {
			std::format_to(out, "[{},{}]", edges[i].first, edges[i].second);
			if (i != edges.size() - 1) {
				std::format_to(out, ",");
			}
		}
	}


	void printPauliCollection(auto out, const auto& collection) {
		std::format_to(out, "    {{\n      \"operators\": [");

		for (size_t i = 0; i < collection.paulis.size(); ++i) {
			std::format_to(out, "\"{}\"", collection.paulis[i]);
			if (i != collection.paulis.size() - 1) {
				std::format_to(out, ",");
			}
		}
		std::format_to(out, "],\n      \"edges\": [");
		printEdgeList(out, collection.graph.getEdges());
		std::format_to(out, "],\n      \"cliffords\": [");

		for (size_t i = 0; i < collection.singleQubitLayer.size(); ++i) {
			const auto& gate = collection.singleQubitLayer[i];
			if (gate == Q::BinaryCliffordGates::I) std::format_to(out, "\"I\"");
			if (gate == Q::BinaryCliffordGates::H) std::format_to(out, "\"H\"");
			if (gate == Q::BinaryCliffordGates::S) std::format_to(out, "\"S\"");
			if (gate == Q::BinaryCliffordGates::SH) std::format_to(out, "\"SH\"");
			if (gate == Q::BinaryCliffordGates::HSH) std::format_to(out, "\"HSH\"");
			if (gate == Q::BinaryCliffordGates::HS) std::format_to(out, "\"HS\"");
			if (i != collection.singleQubitLayer.size() - 1) {
				std::format_to(out, ",");
			}
		}
		std::format_to(out, "]\n    }}");
	}


	void printPauliCollections(auto out, const auto& collections, const MetaInfo& metaInfo) {

		std::format_to(out, "{{\n");
		std::format_to(out, "  \"runtime [seconds]\": {},\n", metaInfo.timeInSeconds);
		std::format_to(out, "  \"num graphs\": {},\n", metaInfo.numGraphs);
		std::format_to(out, "  \"connectivity\": [");
		printEdgeList(out, metaInfo.connectivity.getEdges());
		std::format_to(out, "],\n", metaInfo.numGraphs);
		std::format_to(out, "  \"random seed\": {},\n", metaInfo.randomSeed);
		std::format_to(out, "  \"R_hat_HT\": {},\n", metaInfo.Rhat_HT);
		std::format_to(out, "  \"R_hat_TPB\": {},\n", metaInfo.Rhat_TPB);
		std::format_to(out, "  \"num groups\": {},\n", collections.size());
		std::format_to(out, "  \"num paulis\": {},\n", std::accumulate(collections.begin(), collections.end(), size_t{ 0 }, [](size_t sum, const auto& collection) {return sum+ collection.size(); }));

		//auto mat = metaInfo.connectivity.getAdjacencyMatrix();
		//for(int i=0; i < )

		std::format_to(out, "  \"grouping\": [\n");
		for (size_t i = 0; i < collections.size(); ++i) {
			printPauliCollection(out, collections[i]);
			if (i != collections.size() - 1) {
				std::format_to(out, ",\n");
			}
		}
		std::format_to(out, "\n  ]\n}}\n");
	}


}