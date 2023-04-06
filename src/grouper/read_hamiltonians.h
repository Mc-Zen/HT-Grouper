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
		if (!file) throw std::runtime_error("Error, could not open file");

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
				const auto pauli = Pauli{ trim(pauliAndValue[0], " \"") };

				if (hamiltonian.numQubits == 0) {
					hamiltonian.numQubits = pauli.n;
				}
				if (pauli != Pauli{ hamiltonian.numQubits }) {
					hamiltonian.operators.emplace_back(pauli, std::stod(pauliAndValue[1]));
				}

			}
			hamiltonians.emplace_back(hamiltonian);
		}
		return hamiltonians;
	}
}
