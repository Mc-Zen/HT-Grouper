#pragma once
#include <format>


namespace JsonFormatting {




	void printPauliCollection(auto out, const auto& collection) {
		std::format_to(out, "    {{\n      \"operators\": [");

		for (size_t i = 0; i < collection.paulis.size(); ++i) {
			std::format_to(out, "\"{}\"", collection.paulis[i]);
			if (i != collection.paulis.size() - 1) {
				std::format_to(out, ",");
			}
		}
		std::format_to(out, "],\n      \"edges\": [");
		const auto edges = collection.graph.getEdges();
		for (size_t i = 0; i < edges.size(); ++i) {
			std::format_to(out, "[{},{}]", edges[i].first, edges[i].second);
			if (i != edges.size() - 1) {
				std::format_to(out, ",");
			}
		}
		std::format_to(out, "]\n    }}");
	}


	void printPauliCollections(auto out, const auto& collections) {

		std::format_to(out, "{{\n  \"grouping\": [\n");
		for (size_t i = 0; i < collections.size(); ++i) {
			printPauliCollection(out, collections[i]);
			if (i != collections.size() - 1) {
				std::format_to(out, ",\n");
			}
		}
		std::format_to(out, "\n  ]\n}}\n");
	}


}