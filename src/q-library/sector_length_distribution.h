
#pragma once
#include "graph.h"
#include "special_math.h"
#include "formatting.h"

namespace Q {


	template<int numQubits>
	auto sectorLengthDistribution(const Graph<numQubits>& graph) {
		std::array<int, numQubits + 1> sld{};
		constexpr auto end = pow2(numQubits);

		efficient::BinaryRowMatrix<numQubits, numQubits> adjacencyMatrix{ graph.adjacencyMatrix };

		for (uint64_t i = 0; i < end; ++i) {
			const auto result = adjacencyMatrix * efficient::BinaryVector<numQubits>{i};
			//print("{:0<5b}, {:0<5b} -> {}\n", i, i | result, std::popcount(i | result));
			++sld[(i | result).bitCount()];
		}
		return sld;
	}

}
