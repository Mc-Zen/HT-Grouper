#pragma once
#include <format>


namespace PythonFormatting {




	void printPauliCollection(auto out, const auto& collection) {
		std::format_to(out, "{{ \"operators\": [");
		for (const auto& pauli : collection.paulis) {
			std::format_to(out, "\"{}\",", pauli);
		}
		std::format_to(out, "], \"edges\": [");
		for (const auto& edge : collection.graph.getEdges()) {
			std::format_to(out, "{},", edge);
		}
		std::format_to(out, "]}},\n");
	}
	

	void printPauliCollections(auto out, const auto& collections) {

		std::format_to(out, "grouping = [\n");
		for (const auto& collection : collections) {
			printPauliCollection(out, collection);
			//std::format_to(out, "{} -> {}", collection.paulis, collection.graph.getEdges());
		}
		std::format_to(out, "]\n");
	}


}