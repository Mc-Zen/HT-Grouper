#pragma once
#include "graph.h"
#include "binary_pauli.h"


namespace Q {

	template<int numQubits>
	auto getStabilizer(const Graph<numQubits>& graph) {
		BinaryOperatorSet<numQubits, numQubits> stabilizer;
		for (int row = 0; row < numQubits; ++row) {
			stabilizer[row].z(row) = 1;
			for (int col = 0; col < numQubits; ++col) {
				stabilizer[row].x(col) = graph.adjacencyMatrix(row, col);
			}
		}
		return stabilizer;
	}


	template<int numQubits>
	auto expandStabilizer(const BinaryOperatorSet<numQubits, numQubits>& stabilizer) {
		BinaryOperatorSet<numQubits, pow2(numQubits)> expandedStabilizer{};
		for (int i = 0; i < pow2(numQubits); ++i) {
			for (int j = 0; j < numQubits; ++j) {
				if (i & (1 << j)) {
					expandedStabilizer[i] *= stabilizer[j];
				}
			}
		}
		return expandedStabilizer;
	}

}