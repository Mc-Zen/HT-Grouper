#pragma once

#include <algorithm>
#include "hamiltonian.h"


namespace Q {

	double estimated_shot_reduction(const Hamiltonian& hamiltonian, const std::vector<CollectionWithGraph>& grouping) {
		double numerator{};
		double denominator{};

		for (const auto& group : grouping) {
			double denominatorTerm{};
			for (const auto& pauli : group.paulis) {
				auto coefficient = std::get<double>(*std::ranges::find(hamiltonian.operators, pauli, [](const auto& a) { return std::get<0>(a); }));
				double absolute = std::abs(coefficient);
				numerator += absolute;
				denominatorTerm += absolute * absolute;
			}
			denominator += std::sqrt(denominatorTerm);
		}
		return numerator * numerator / (denominator * denominator);
	}
}
