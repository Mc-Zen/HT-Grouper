#pragma once
#include <fstream>
#include <string>
#include "string_utility.h"

namespace Q {


	class ConfigReadError : public std::runtime_error {
	public:
		using std::runtime_error::runtime_error;
	};

	struct Configuration {
		std::string filename;
		std::string outfilename;
		std::string connectivity;
		int numThreads{};
		int maxEdgeCount{};
		int numGraphs{};
		bool sortGraphsByEdgeCount{};
	};

	Configuration readConfig(const std::string& filename) {

		std::ifstream file{ filename };
		if (!file) throw ConfigReadError(std::format("Could not open file \"{}\"", filename));

		Configuration config;

		std::string line;
		while (std::getline(file, line)) {
			if (line.empty()) continue;
			line = trim(split(line, '#')[0], " \t"); // strip comments and whitespace
			if (line.empty()) continue;
			auto components = split(line, '=');
			if (components.size() != 2) throw ConfigReadError(std::format("Invalid attribute format for attribute \"{}\". Name and value need to be seperated by a \"=\" sign. ", line));

			auto name = trim(components[0], " \t");
			auto value = trim(components[1], " \t");

			if (name == "filename") {
				if (config.filename != "") throw ConfigReadError("Duplicate attribute \"filename\"");
				config.filename = value;
			}
			else if (name == "outfilename") {
				if (config.outfilename != "") throw ConfigReadError("Duplicate attribute \"outfilename\"");
				config.outfilename = value;
			}
			else if (name == "connectivity") {
				if (config.connectivity != "") throw ConfigReadError("Duplicate attribute \"connectivity\"");
				config.connectivity = value;
			}
			else if (name == "numThreads") {
				if (config.numThreads != 0) throw ConfigReadError("Duplicate attribute \"numThreads\"");
				int numThreads = std::stoi(value);
				if (numThreads < 1 || numThreads > 255) throw ConfigReadError("The \"numThreads\" attribute can only take values between 1 and 255");
				config.numThreads = numThreads;
			}
			else if (name == "maxEdgeCount") {
				if (config.maxEdgeCount != 0) throw ConfigReadError("Duplicate attribute \"maxEdgeCount\"");
				int maxEdgeCount = std::stoi(value);
				if (maxEdgeCount < 1) throw ConfigReadError("The \"maxEdgeCount\" attribute needs to be positive");
				config.maxEdgeCount = maxEdgeCount;
			}

			else if (name == "numGraphs") {
				if (config.numGraphs != 0) throw ConfigReadError("Duplicate attribute \"numGraphs\"");
				int numGraphs = std::stoi(value);
				if (numGraphs < 1) throw ConfigReadError("The \"numGraphs\" attribute needs to be positive");
				config.numGraphs = numGraphs;
			}
			else if (name == "sortGraphsByEdgeCount") {
				bool sortGraphsByEdgeCount;
				if (value == "true") sortGraphsByEdgeCount = true;
				else if (value == "false") sortGraphsByEdgeCount = false;
				else throw ConfigReadError("The \"sortGraphsByEdgeCount\" attribute can only be true or false");
				config.sortGraphsByEdgeCount = sortGraphsByEdgeCount;
			}
			else {
				throw ConfigReadError(std::format("Unknown attribute \"{}\"", name));
			}
		}
		return config;
	}


	class ConnectivityError : public std::runtime_error {
	public:
		using std::runtime_error::runtime_error;
	};



	class Connectivity {
	public:
		enum class Type { Linear, Cycle, Star, Matrix };
		using AdjacencyMatrix = Math::Matrix<Q::Binary>;

		Connectivity(Type type) : type(type) {}
		Connectivity(AdjacencyMatrix matrix) : adjacencyMatrix(matrix), type(Type::Matrix) {}

		Graph<> getGraph(int numQubits) const {
			switch (type) {
				using enum Type;
			case Linear: return Graph<>::linear(numQubits);
			case Cycle: return Graph<>::cycle(numQubits);
			case Star: return Graph<>::star(numQubits);
			}
			if (numQubits != adjacencyMatrix.rows()) throw ConnectivityError(std::format("The adjacency matrix has {} qubits while {} were specified", adjacencyMatrix.rows(), numQubits));
			auto graph = Graph<>(numQubits);
			graph.adjacencyMatrix = adjacencyMatrix;
			return graph;

		}

	private:
		Type type;
		AdjacencyMatrix adjacencyMatrix;
	};

	Connectivity readConnectivity(const std::string& filename) {

		std::ifstream file{ filename };
		if (!file) throw ConnectivityError(std::format("Could not open file \"{}\"", filename));


		std::string line;
		Math::Matrix<Q::Binary> matrix;
		int matrixLine{};
		while (std::getline(file, line)) {
			if (line.empty()) continue;
			line = trim(split(line, '#')[0], " \t"); // strip comments and whitespace
			if (line.empty()) continue;
			if (line == "linear")
				return Connectivity{ Connectivity::Type::Linear };
			if (line == "cycle")
				return Connectivity{ Connectivity::Type::Cycle };
			if (line == "star")
				return Connectivity{ Connectivity::Type::Star };

			auto entries = split(line, ' ');
			int numQubits = entries.size();
			if (matrix.rows() == 0) {
				matrix.resize(numQubits, numQubits);
			}
			else if (matrix.rows() != numQubits) {
				throw ConnectivityError(std::format("Each row of the matrix needs to have the same number of entries, row \"{}\"", line));
			}

			for (int i = 0; i < numQubits; ++i) {
				matrix(matrixLine, i) = std::stoi(entries[i]);
			}

			++matrixLine;

		}
		return Connectivity{ matrix };
	}

}