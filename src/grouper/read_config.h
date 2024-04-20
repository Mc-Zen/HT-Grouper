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
		int64_t numThreads{};
		int64_t maxEdgeCount{};
		int64_t numGraphs{};
		int64_t intermediateFileFrequency{};
		bool sortGraphsByEdgeCount{ true };
		bool extractComputationalBasis{ true };
		unsigned int seed{};
	};


	int64_t string_to_int(const std::string& str) {
		try {
			return std::stoll(str);
		}
		catch (std::out_of_range& e) {
			throw ConfigReadError(fmt::format("Integer out of range: \"{}\"", str));
		}
	}

	Configuration readConfig(const std::string& filename) {

		std::ifstream file{ filename };
		if (!file) throw ConfigReadError(fmt::format("Could not open file \"{}\"", filename));

		Configuration config;

		std::string line;
		while (std::getline(file, line)) {
			if (line.empty()) continue;
			line = trim(split(line, '#')[0], " \t"); // strip comments and whitespace
			if (line.empty()) continue;
			auto components = splitOnce(line, '=');
			if (components.size() != 2) throw ConfigReadError(fmt::format("Invalid attribute format for attribute \"{}\". Name and value need to be seperated by a \"=\" sign. ", line));

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
				auto numThreads = string_to_int(value);
				if (numThreads < 1 || numThreads > 255) throw ConfigReadError("The \"numThreads\" attribute can only take values between 1 and 255");
				config.numThreads = numThreads;
			}
			else if (name == "maxEdgeCount") {
				if (config.maxEdgeCount != 0) throw ConfigReadError("Duplicate attribute \"maxEdgeCount\"");
				auto maxEdgeCount = string_to_int(value);
				if (maxEdgeCount < 1) throw ConfigReadError("The \"maxEdgeCount\" attribute needs to be positive");
				config.maxEdgeCount = maxEdgeCount;
			}

			else if (name == "numGraphs") {
				if (config.numGraphs != 0) throw ConfigReadError("Duplicate attribute \"numGraphs\"");
				auto numGraphs = string_to_int(value);
				if (numGraphs < 1) throw ConfigReadError("The \"numGraphs\" attribute needs to be positive");
				config.numGraphs = numGraphs;
			}
			else if (name == "seed") {
				if (config.seed != 0) throw ConfigReadError("Duplicate attribute \"seed\"");
				auto seed = string_to_int(value);
				if (seed < 1) throw ConfigReadError("The \"seed\" attribute needs to be positive");
				config.seed = seed;
			}
			else if (name == "sortGraphsByEdgeCount") {
				bool sortGraphsByEdgeCount;
				if (value == "true") sortGraphsByEdgeCount = true;
				else if (value == "false") sortGraphsByEdgeCount = false;
				else throw ConfigReadError("The \"sortGraphsByEdgeCount\" attribute can only be true or false");
				config.sortGraphsByEdgeCount = sortGraphsByEdgeCount;
			}
			else if (name == "extractComputationalBasis") {
				bool extractComputationalBasis;
				if (value == "true") extractComputationalBasis = true;
				else if (value == "false") extractComputationalBasis = false;
				else throw ConfigReadError("The \"extractComputationalBasis\" attribute can only be true or false");
				config.extractComputationalBasis = extractComputationalBasis;
			}
			else if (name == "intermediateFileFrequency") {
				if (config.intermediateFileFrequency != 0) throw ConfigReadError("Duplicate attribute \"intermediateFileFrequency\"");
				auto intermediateFileFrequency = string_to_int(value);
				if (intermediateFileFrequency < 0) throw ConfigReadError("The \"intermediateFileFrequency\" attribute needs to be positive or zero");
				config.intermediateFileFrequency = intermediateFileFrequency;
			}
			else {
				throw ConfigReadError(fmt::format("Unknown attribute \"{}\"", name));
			}
		}

		if (config.filename == "")
			throw ConfigReadError("No [filename] specified");
		if (config.outfilename == "")
			throw ConfigReadError("No [outfilename] specified");
		if (config.connectivity == "")
			throw ConfigReadError("No [connectivity] specified");
		if (config.numGraphs == 0) config.numGraphs = 100;
		if (config.maxEdgeCount == 0) config.maxEdgeCount = 1000;
		if (config.numThreads == 0) config.numThreads = 1;

		return config;
	}


	class ConnectivityError : public std::runtime_error {
	public:
		using std::runtime_error::runtime_error;
	};




	class Connectivity {
	public:
		enum class Type { Linear, Cycle, Star, Matrix, All, SquareLattice };
		using AdjacencyMatrix = Math::Matrix<Q::Binary>;

		Connectivity(Type type) : type(type) {}
		Connectivity(AdjacencyMatrix matrix) : adjacencyMatrix(matrix), type(Type::Matrix) {}

		Graph<> getGraph(int numQubits) const {
			switch (type) {
				using enum Type;
			case Linear: return Graph<>::linear(numQubits);
			case Cycle: return Graph<>::cycle(numQubits);
			case Star: return Graph<>::star(numQubits);
			case All: return Graph<>::fullyConnected(numQubits);
			case SquareLattice: return Graph<>::squareLattice(numQubits);
			}
			if (numQubits != adjacencyMatrix.rows()) throw ConnectivityError(fmt::format("The adjacency matrix has {} qubits while {} were specified", adjacencyMatrix.rows(), numQubits));
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
		if (!file) throw ConnectivityError(fmt::format("Could not open file \"{}\"", filename));


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
			if (line == "all")
				return Connectivity{ Connectivity::Type::All };
			if (line == "square-lattice")
				return Connectivity{ Connectivity::Type::SquareLattice };

			auto entries = split(line, ' ');
			int numQubits = entries.size();
			if (matrix.rows() == 0) {
				matrix.resize(numQubits, numQubits);
			}
			else if (matrix.rows() != numQubits) {
				throw ConnectivityError(fmt::format("Each row of the matrix needs to have the same number of entries, row \"{}\"", line));
			}

			for (int i = 0; i < numQubits; ++i) {
				matrix(matrixLine, i) = std::stoi(entries[i]);
			}

			++matrixLine;

		}
		return Connectivity{ matrix };
	}


}