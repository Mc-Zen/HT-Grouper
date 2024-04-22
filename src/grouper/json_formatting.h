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
			fmt::format_to(out, "[{},{}]", edges[i].first, edges[i].second);
			if (i != edges.size() - 1) {
				fmt::format_to(out, ",");
			}
		}
	}


	void printPauliCollection(auto out, const auto& collection) {
		fmt::format_to(out, "    {{\n      \"operators\": [");

		for (size_t i = 0; i < collection.paulis.size(); ++i) {
			fmt::format_to(out, "\"{}\"", collection.paulis[i]);
			if (i != collection.paulis.size() - 1) {
				fmt::format_to(out, ",");
			}
		}
		fmt::format_to(out, "],\n      \"edges\": [");
		printEdgeList(out, collection.graph.getEdges());
		fmt::format_to(out, "],\n      \"cliffords\": [");

		for (size_t i = 0; i < collection.singleQubitLayer.size(); ++i) {
			const auto& gate = collection.singleQubitLayer[i];
			if (gate == Q::BinaryCliffordGates::I) fmt::format_to(out, "\"I\"");
			if (gate == Q::BinaryCliffordGates::H) fmt::format_to(out, "\"H\"");
			if (gate == Q::BinaryCliffordGates::S) fmt::format_to(out, "\"S\"");
			if (gate == Q::BinaryCliffordGates::SH) fmt::format_to(out, "\"SH\"");
			if (gate == Q::BinaryCliffordGates::HSH) fmt::format_to(out, "\"HSH\"");
			if (gate == Q::BinaryCliffordGates::HS) fmt::format_to(out, "\"HS\"");
			if (i != collection.singleQubitLayer.size() - 1) {
				fmt::format_to(out, ",");
			}
		}
		fmt::format_to(out, "]\n    }}");
	}


	void printPauliCollections(auto out, const auto& collections, const MetaInfo& metaInfo) {

		size_t numPaulis = std::accumulate(collections.begin(), collections.end(), size_t{ 0 }, [](size_t c, const auto& col) { return c + col.size(); });
		fmt::format_to(out, "{{\n");
		fmt::format_to(out, "  \"runtime [seconds]\": {},\n", metaInfo.timeInSeconds);
		fmt::format_to(out, "  \"num graphs\": {},\n", metaInfo.numGraphs);
		fmt::format_to(out, "  \"num paulis\": {},\n", numPaulis);
		fmt::format_to(out, "  \"num groups\": {},\n", collections.size());
		fmt::format_to(out, "  \"connectivity\": [");
		printEdgeList(out, metaInfo.connectivity.getEdges());
		fmt::format_to(out, "],\n", metaInfo.numGraphs);
		fmt::format_to(out, "  \"random seed\": {},\n", metaInfo.randomSeed);
		fmt::format_to(out, "  \"R_hat_HT\": {},\n", metaInfo.Rhat_HT);
		fmt::format_to(out, "  \"R_hat_TPB\": {},\n", metaInfo.Rhat_TPB);

		//auto mat = metaInfo.connectivity.getAdjacencyMatrix();
		//for(int i=0; i < )

		fmt::format_to(out, "  \"grouping\": [\n");
		for (size_t i = 0; i < collections.size(); ++i) {
			printPauliCollection(out, collections[i]);
			if (i != collections.size() - 1) {
				fmt::format_to(out, ",\n");
			}
		}
		fmt::format_to(out, "\n  ]\n}}\n");
	}


}