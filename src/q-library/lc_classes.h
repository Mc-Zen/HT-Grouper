#pragma once

#include "binary_pauli.h"
#include "mub.h"
#include "graph.h"
#include "efficient_mub.h"
#include "stabilizer.h"
#include <vector>
#include <variant>
#include <format>

namespace Q {

	struct Single {
		int v{};
	};

	struct Pair {
		Pair(int v1, int v2) : v1(std::min(v1, v2)), v2(std::max(v1, v2)) {}
		int v1{};
		int v2{};
		constexpr friend bool operator==(const Pair& a, const Pair& b) = default;
	};

	struct Triple {
		Triple(int v1, int v2, int v3) : v1(std::min({ v1, v2, v3 })), v3(std::max({ v1, v2, v3 })) {
			this->v2 = v1 + v2 + v3 - this->v1 - this->v3;
		}
		int v1{};
		int v2{};
		int v3{};
		constexpr friend bool operator==(const Triple& a, const Triple& b) = default;

	};

	struct PairAndTriple {
		Pair pair;
		Triple triple;
		constexpr friend bool operator==(const PairAndTriple& a, const PairAndTriple& b) = default;
	};


	struct TwoPairs {
		Pair pair1;
		Pair pair2;
		constexpr friend bool operator==(const TwoPairs& a, const TwoPairs& b) = default;
	};

	template<int n>
	struct LCClass;

	template<>
	struct LCClass<2> {
		enum class Type {
			Separable,  // Separable graph
			Entangled,  // The vertices are connected by an edge
		};

		Type type;
	};

	template<>
	struct LCClass<3> {
		enum class Type {
			Separable,  // Separable graph
			Pair,		// The vertices are connected by an edge
			Triple,     // One triple of connected vertices, the remaining vertex is isolated
		};

		Type type;
		std::variant<std::monostate, Pair> pair;
	};

	template<>
	struct LCClass<4> {
		enum class Type {
			Separable,  // Separable graph
			Pair,		// The vertices are connected by an edge
			Triple,     // One triple of connected vertices, the remaining vertex is isolated
			TwoPairs,	// Two pairs of connected vertices with no edges between them, the remaining vertex is isolated
			Star,		// A star of 4 vertices
			Line,		// A path graph of 4 vertices
		};

		Type type;
		std::variant<std::monostate, Pair, Triple, TwoPairs> data;
	};

	template<>
	struct LCClass<5> {
		enum class Type {
			Separable,		// Separable graph
			Pair,			// The vertices are connected by an edge
			Triple,			// One triple of connected vertices, the remaining vertex is isolated
			TwoPairs,		// Two pairs of connected vertices with no edges between them, the remaining vertex is isolated
			Star4,			// A star of 4 vertices
			Line4,			// A path graph of 4 vertices
			Star,			// A star of 5 vertices
			PairAndTriple,	// One pair and one triple of connected vertices with no edges between them
			T,				// A star of four vertices where the remaining vertex is connected to one of the former
			Line,			// A path graph of 5 vertices
			Cycle,			// A cycle graph
		};

		Type type;
		std::variant<std::monostate, Single, Pair, Triple, TwoPairs, PairAndTriple> data;
	};



	template<int n>
	LCClass<n> determine_lc_class(const BinaryOperatorSet<n, n>& stabilizer);

	template<>
	LCClass<2> determine_lc_class(const BinaryOperatorSet<2, 2>& stabilizer) {
		auto areQubitsEntangled = Q::areQubitsEntangled(stabilizer);
		auto numEntangledQubits = std::ranges::count(areQubitsEntangled, true);
		if (numEntangledQubits == 0) {
			return LCClass<2>{LCClass<2>::Type::Separable};
		}
		else if (numEntangledQubits == 2) {
			return LCClass<2>{LCClass<2>::Type::Entangled};
		}
		throw std::runtime_error("Invalid stabilizer");
	}

	template<int n>
	auto getEntanglementInformation(const BinaryOperatorSet<n, n>& stabilizer) {
		const auto areQubitsEntangled = Q::areQubitsEntangled(stabilizer);
		const auto numEntangledQubits = std::ranges::count(areQubitsEntangled, true);
		std::vector<int> entangledQubits;
		for (int i = 0; i < n; ++i) {
			if (areQubitsEntangled[i]) entangledQubits.push_back(i);
		}
		return std::make_tuple(areQubitsEntangled, numEntangledQubits, entangledQubits);
	}

	template<>
	LCClass<3> determine_lc_class(const BinaryOperatorSet<3, 3>& stabilizer) {
		constexpr int n = 3;
		const auto [areQubitsEntangled, numEntangledQubits, entangledQubits] = getEntanglementInformation(stabilizer);

		switch (numEntangledQubits) {
		case 0:
			return LCClass<n>{LCClass<n>::Type::Separable};
		case 2:
			return LCClass<n>{LCClass<n>::Type::Pair, Pair{ entangledQubits[0],entangledQubits[1] }};
		case 3:
			return LCClass<n>{LCClass<n>::Type::Triple};
		default:
			throw std::runtime_error("Invalid stabilizer");
		}
	}

	template<>
	LCClass<4> determine_lc_class(const BinaryOperatorSet<4, 4>& stabilizer) {
		constexpr int n = 4;
		const auto [areQubitsEntangled, numEntangledQubits, entangledQubits] = getEntanglementInformation(stabilizer);

		switch (numEntangledQubits) {
		case 0:
			return LCClass<n>{LCClass<n>::Type::Separable};
		case 2:
			return LCClass<n>{LCClass<n>::Type::Pair, Pair{ entangledQubits[0],entangledQubits[1] }};
		case 3:
			return LCClass<n>{LCClass<n>::Type::Triple, Triple{ entangledQubits[0], entangledQubits[1],entangledQubits[2] } };
		case 4:
		{
			//const auto fullStabilizer = constructMubSetFromCanonicalGeneratingSet(stabilizer);
			const auto fullStabilizer = efficient::expandStabilizer(efficient::toEfficientStabilizer(stabilizer));
			bool has_IIAA = efficient::countIdentityStructure(fullStabilizer, { 0b1100 }) != 0;
			bool has_IAIA = efficient::countIdentityStructure(fullStabilizer, { 0b1010 }) != 0;
			bool has_AIIA = efficient::countIdentityStructure(fullStabilizer, { 0b0110 }) != 0;
			if (has_IIAA && has_AIIA) {
				return LCClass<n>{LCClass<n>::Type::Star};
			}

			const int A_n = efficient::countIdentityStructure(fullStabilizer, { 0b0000 }); // last entry of the sector length distribution


			TwoPairs twoPairs{ {0,0},{0,0} };
			if (has_IIAA) {
				twoPairs = { {0,1},{2,3} };
			}
			else if (has_IAIA) {
				twoPairs = { {0,2},{1,3} };
			}
			else if (has_AIIA) {
				twoPairs = { {0,3},{1,2} };
			}
			else {
				throw std::runtime_error("Invalid stabilizer");
			}
			if (A_n == 9) {
				return LCClass<n>{LCClass<n>::Type::TwoPairs, twoPairs};
			}
			else if (A_n == 5) {
				return LCClass<n>{LCClass<n>::Type::Line, twoPairs};
			}
		}
		default: break;
		}
		throw std::runtime_error("Invalid stabilizer");
	}

}
