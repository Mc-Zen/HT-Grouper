#pragma once

#include <algorithm>
#include "hamiltonian.h"


namespace Q {

	/// @brief Compute estimated shot reduction compared to single Pauli measurements. 
	/// 
	///           ⎛  ∑_i^N ∑_j^{m_i} |a_ij|    ⎞²
	/// \hat{R} = ⎜----------------------------⎟  
	///	          ⎝ ∑_i^N √(∑_j^{m_i} |a_ij|²) ⎠
	/// 
	/// as defined in https://doi.org/10.22331/q-2021-01-20-385
	/// 
	/// @param hamiltonian Hamiltonian 
	/// @param grouping    Grouping of the operators in hamiltonian
	/// @return            Estimated shot reduction
	double estimated_shot_reduction(const Hamiltonian& hamiltonian, const std::vector<CollectionWithGraph>& grouping) {
		double numerator{};
		double denominator{};

		for (const auto& group : grouping) {
			double denominatorTerm{};
			for (const auto& pauli : group.paulis) {
				if (pauli == Pauli::Identity(hamiltonian.numQubits)) continue; // no need to measure identity

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
