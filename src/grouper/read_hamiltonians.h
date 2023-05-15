#pragma once

#include <fstream>
#include <string>
#include "formatting.h"
#include "binary_pauli.h"
#include "string_utility.h"
#include "hamiltonian.h"

namespace Q {


	/// @brief Read hamiltonians from python file in form of a dictionary
	/// @param filename Path to file
	/// @return List of hamiltonian specifications
	std::vector<Hamiltonian> readHamiltonians(const std::string& filename) {

		std::ifstream file{ filename };
		if (!file) throw std::runtime_error(std::format("Error, could not open file {}", filename));

		std::vector<Hamiltonian> hamiltonians;

		std::string line;
		while (std::getline(file, line)) {
			if (line.empty()) continue;
			auto dictStart = line.find('{');
			auto dictEnd = line.rfind('}');
			if (dictStart == std::string::npos || dictEnd == std::string::npos) throw std::runtime_error("Error, wrong format");

			line = line.substr(dictStart + 1, dictEnd - dictStart - 1);

			auto paulisStart = line.find('{');
			auto paulisEnd = line.rfind('}');
			if (paulisStart == std::string::npos || paulisEnd == std::string::npos) throw std::runtime_error("Error, wrong format");
			auto pauliStrings = split(line.substr(paulisStart + 1, paulisEnd - paulisStart - 1), ',');

			Hamiltonian hamiltonian;
			for (const auto& pauliString : pauliStrings) {
				auto pauliAndValue = split(pauliString, ':');
				if (pauliAndValue.size() != 2) throw std::runtime_error("Error, wrong format");
				const auto pauli = Pauli{ trim(pauliAndValue[0], " \"\'") };

				if (hamiltonian.numQubits == 0) {
					hamiltonian.numQubits = pauli.numQubits();
				}
				if (pauli != Pauli{ hamiltonian.numQubits }) {
					hamiltonian.operators.emplace_back(pauli, std::stod(pauliAndValue[1]));
				}

			}
			hamiltonians.emplace_back(hamiltonian);
		}
		return hamiltonians;
	}





	class ReadHamiltonianError : public std::runtime_error {
	public:
		using std::runtime_error::runtime_error;
	};

	/// @brief Read hamiltonians from python file in form of a dictionary
	/// @param filename Path to file
	/// @return List of hamiltonian specifications
	Hamiltonian readHamiltonianFromJson(const std::string& filename) {

		std::ifstream file{ filename };
		if (!file) throw ReadHamiltonianError(std::format("Error, could not open file {}", filename));

		Hamiltonian hamiltonian;

		std::string line;
		int lineIndex{ 0 };
		while (std::getline(file, line)) {
			++lineIndex;
			if (line.empty()) continue;
			line = trim(line, " \t{}");
			if (line.empty()) continue;

			auto components = split(line, ':');
			if (components.size() != 2) throw ReadHamiltonianError(std::format("Invalid format for \"{}\" at line {}", line, lineIndex));

			auto pauliString = trim(components[0], " \t\"\'");
			auto value = trim(components[1], " \t,");

			if (pauliString.size() == 0) throw ReadHamiltonianError(std::format("Empty Pauli string at line {}", lineIndex));
			Pauli pauli{ pauliString };
			if (hamiltonian.numQubits == 0) {
				hamiltonian.numQubits = pauli.numQubits();
			}
			else if (hamiltonian.numQubits != pauli.numQubits()) {
				throw ReadHamiltonianError(std::format("The Pauli {} at line {} does not have the same number of qubits as the preceding Paulis", pauliString, lineIndex));
			}

			try {
				auto coefficient = std::stod(value);
				hamiltonian.operators.emplace_back(pauli, coefficient);
			}
			catch (std::invalid_argument&) {
				throw ReadHamiltonianError(std::format("Invalid coefficient {} at line {}", value, lineIndex));
			}
			catch (std::out_of_range&) {
				throw ReadHamiltonianError(std::format("Out of range coefficient {} at line {}", value, lineIndex));
			}
		}
		return hamiltonian;
	}



	/// @brief Read Pauli groups from file, in the following format:
	///        {XYZ,ZZX,IXX}
	///        {XZZ}
	///        {XZX,YZY,YYY,IIY,IXI}
	///        ...
	/// @param filename Path to file
	/// @return List of Pauli groups
	std::vector<std::vector<Pauli>> readPauliGroups(const std::string& filename) {

		std::ifstream file{ filename };
		if (!file) throw std::runtime_error(std::format("Error, could not open file {}", filename));

		std::vector<std::vector<Pauli>> groups;

		std::string line;
		while (std::getline(file, line)) {
			if (line.empty()) continue;

			std::vector<Pauli> group;
			line = trim(line, "{} ");
			auto paulis = split(line, ",");
			for (const auto& pauli : paulis) {
				group.emplace_back(pauli);
			}
			groups.emplace_back(std::move(group));
		}
		return groups;
	}
}
