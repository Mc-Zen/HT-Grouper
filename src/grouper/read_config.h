#pragma once
#include <fstream>
#include <string>
#include <charconv>
#include <variant>
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
		bool generateTPBs{ true };
		bool extractComputationalBasis{ true };
		unsigned int seed{};
		bool verboseLog{ true };
		int64_t grouperType{ 1 };
	};

	struct IntValue {
		IntValue() = default;
		IntValue(int64_t initial, int64_t min, int64_t max) : value(initial), min(min), max(max) {}
		void read(std::string_view input) {
			auto result = std::from_chars(input.data(), input.data() + input.size(), value);
			if (result.ec == std::errc::invalid_argument) { throw ConfigReadError(fmt::format("Integer out of range: \"{}\"", input)); }
			if (value < min || value > max) { throw ConfigReadError(fmt::format("Integer needs to be in range [{}, {}] (was {})", min, max, value)); }
		}

		auto write() const { return std::to_string(value); }

		int64_t value{ 0 };
		int64_t min{ std::numeric_limits<int64_t>::min() };
		int64_t max{ std::numeric_limits<int64_t>::max() };
	};

	struct StringValue {
		explicit StringValue(std::string_view initial = "") : value(initial) {};
		void read(std::string_view input) { value = input; }
		std::string write() const { return '"' + value + '"'; }

		std::string value;
	};

	struct BoolValue {
		explicit BoolValue(bool initial = false) : value(initial) {};
		void read(std::string_view input) {
			if (toLower(input) == "true") { value = true; }
			else if (toLower(input) == "false") { value = false; }
			else { throw ConfigReadError(fmt::format("Expected bool value, got {}", input)); }
		}
		std::string write() const { return value ? "true" : "false"; }
		bool value;
	};


	class Attribute {
	public:

		template<class K>
		Attribute(std::string_view name, K value) : name_(name), value_(value) {}

		void read(std::string_view input) {
			if (std::holds_alternative<IntValue>(value_)) { value<IntValue>().read(input); }
			else if (std::holds_alternative<BoolValue>(value_)) { value<BoolValue>().read(input); }
			else if (std::holds_alternative<StringValue>(value_)) { value<StringValue>().read(input); }
		}

		std::string write() const {
			if (std::holds_alternative<IntValue>(value_)) { return value<IntValue>().write(); }
			if (std::holds_alternative<BoolValue>(value_)) { return value<BoolValue>().write(); }
			if (std::holds_alternative<StringValue>(value_)) { return value<StringValue>().write(); }
			return {};
		}


		std::string_view name() const { return name_; }

		template<class T>
		const T& value() const { return std::get<T>(value_); }
		template<class T>
		T& value() { return std::get<T>(value_); }

	private:
		std::string name_;
		std::variant<IntValue, StringValue, BoolValue> value_;
	};


	//template<class T>
	//struct MakeAttributeBase {
	//	std::string_view name;
	//	T initial{ 0 };
	//};

	template<class T>
	struct MakeAttribute {
	};

	template<>
	struct MakeAttribute<int64_t> {
		operator Attribute() const { return Attribute(this->name, IntValue{ this->initial, min, max }); }
		std::string_view name;
		int64_t initial{ 0 };
		int64_t min{ std::numeric_limits<int64_t>::min() };
		int64_t max{ std::numeric_limits<int64_t>::max() };
	};

	template<>
	struct MakeAttribute<bool> {
		operator Attribute() const { return Attribute(this->name, BoolValue{this->initial}); }
		std::string_view name;
		bool initial{};
	};

	template<>
	struct MakeAttribute<std::string> {
		operator Attribute() const { return Attribute(this->name, StringValue{this->initial}); }
		std::string_view name;
		std::string_view initial{};
	};



	struct Config {

		void readAttribute(std::string_view name, std::string_view value) {
			getAttribute(name).read(value);
		}

		template<class T>
		T get(std::string_view name);

		Attribute& getAttribute(std::string_view name) {
			auto attr = std::find_if(attributes.begin(), attributes.end(), [&name](const auto& attribute) { return attribute.name() == name; });
			if (attr == attributes.end()) { throw ConfigReadError(fmt::format("Unknown attribute '{}'", name)); }
			return *attr;
		}

		const Attribute& getAttribute(std::string_view name) const {
			auto attr = std::find_if(attributes.begin(), attributes.end(), [&name](const auto& attribute) { return attribute.name() == name; });
			if (attr == attributes.end()) { throw ConfigReadError(fmt::format("Unknown attribute '{}'", name)); }
			return *attr;
		}


		std::vector<Attribute> attributes{
			MakeAttribute<std::string>{.name = "config" },
			MakeAttribute<std::string>{.name = "filename" },
			MakeAttribute<std::string>{.name = "outfilename" },
			MakeAttribute<std::string>{.name = "connectivity" },
			MakeAttribute<int64_t>{.name = "numThreads", .initial = 1, .min = 1, .max = 10000},
			MakeAttribute<int64_t>{.name = "numGraphs", .initial = 100, .min = 1 },
			MakeAttribute<int64_t>{.name = "maxEdgeCount", .initial = 1000, .min = 0,},
			MakeAttribute<int64_t>{.name = "intermediateFileFrequency", .initial = 0, .min = 0 },
			MakeAttribute<int64_t>{.name = "grouperType", .initial = 1, .min = 1, .max = 2 },
			MakeAttribute<int64_t>{.name = "seed", .initial = 0 },
			MakeAttribute<bool>{.name = "sortGraphsByEdgeCount", .initial = true },
			MakeAttribute<bool>{.name = "generateTPBs", .initial = true },
			MakeAttribute<bool>{.name = "extractComputationalBasis", .initial = true },
			MakeAttribute<bool>{.name = "verboseLog", .initial = true },
		};
	};

	template<>
	int64_t Config::get<int64_t>(std::string_view name) { return getAttribute(name).value<IntValue>().value; }

	template<>
	std::string Config::get<std::string>(std::string_view name) { return getAttribute(name).value<StringValue>().value; }

	template<>
	bool Config::get<bool>(std::string_view name) { return getAttribute(name).value<BoolValue>().value; }


	void fillConfigFromFile(const std::string& filename, Config& config) {
		std::ifstream file{ filename };
		if (!file) throw ConfigReadError(fmt::format("Could not open file \"{}\"", filename));

		std::string line;
		while (std::getline(file, line)) {
			if (line.empty()) continue;
			line = trim(split(line, '#')[0], " \t"); // strip comments and whitespace
			if (line.empty()) continue;
			auto components = splitOnce(line, '=');
			if (components.size() != 2) throw ConfigReadError(fmt::format("Invalid attribute format for attribute \"{}\". Name and value need to be seperated by a \"=\" sign. ", line));

			auto name = trim(components[0], " \t");
			auto value = trim(components[1], " \t");
			config.readAttribute(name, value);
		}
	}

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
			else if (name == "grouperType") {
				if (config.grouperType != 1) throw ConfigReadError("Duplicate attribute \"grouperType\"");
				auto grouperType = string_to_int(value);
				if (grouperType < 1 || grouperType > 2) throw ConfigReadError("The \"grouperType\" attribute can only take values 1 and 2");
				config.grouperType = grouperType;
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
			else if (name == "generateTPBs") {
				bool generateTPBs;
				if (value == "true") generateTPBs = true;
				else if (value == "false") generateTPBs = false;
				else throw ConfigReadError("The \"generateTPBs\" attribute can only be true or false");
				config.generateTPBs = generateTPBs;
			}
			else if (name == "verboseLog") {
				bool verboseLog;
				if (value == "true") verboseLog = true;
				else if (value == "false") verboseLog = false;
				else throw ConfigReadError("The \"verboseLog\" attribute can only be true or false");
				config.verboseLog = verboseLog;
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
		enum class Type { Linear, Cycle, Star, Matrix, All, SquareLattice, Empty };
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
			case Empty: return Graph<>(numQubits);
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
			if (line == "empty")
				return Connectivity{ Connectivity::Type::Empty };

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