#pragma once

#include "pauli.h"
#include <vector>
#include <utility>

namespace Q {

	struct Hamiltonian {
		std::vector<std::pair<Pauli, double>> operators;
		int numQubits{};
	};

}
