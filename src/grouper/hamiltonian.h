#pragma once

#include "pauli.h"

namespace Q {

	struct Hamiltonian {
		std::vector<std::pair<Pauli, double>> operators;
		int numQubits{};
	};

}