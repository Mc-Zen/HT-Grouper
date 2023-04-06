
#pragma once
#include "efficient_binary_math.h"
#include "special_math.h"
#include "random_clifford_element.h"
#include <random>

namespace Q::efficient {


	template<int numQubits>
	struct BinaryPauliOperator {
	public:
		BinaryVector<numQubits> r;
		BinaryVector<numQubits> s;

		constexpr BinaryPauliOperator() = default;
		explicit constexpr BinaryPauliOperator(const Q::BinaryPauliOperator<numQubits>& op) {
			for (size_t i = 0; i < numQubits; ++i) {
				r.set(i, op[i][0]);
				s.set(i, op[i][1]);
			}
		}

		constexpr Q::BinaryPauliOperator<numQubits> toBinaryPauliOperator() const {
			Q::BinaryPauliOperator<numQubits> op;
			for (size_t i = 0; i < numQubits; ++i) {
				op[i][0] = r.get(i);
				op[i][1] = s.get(i);
			}
			return op;
		}

		constexpr BinaryPauliOperator applyClifford(const Clifford<numQubits>& cliff) const {
			BinaryPauliOperator result;
			result.r = cliff.Axx * r + cliff.Axz * s;
			result.s = cliff.Azx * r + cliff.Azz * s;
			return result;
		}

		int pauliWeight() const { return (r | s).bitCount(); }
		int identityCount() const { return numQubits - pauliWeight(); }
		BinaryVector<numQubits> getIdentityString() const { return ~(r | s); }

		std::string toString() const {
			std::string str;
			str.reserve(numQubits);
			for (size_t i = 0; i < numQubits; ++i) str += Q::toChar({ r.get(i), s.get(i) });
			return str;
		}

		constexpr friend bool operator==(const BinaryPauliOperator& a, const BinaryPauliOperator& b) = default;
	};




	template<int numQubits, int m>
	using OperatorSet = std::array<BinaryPauliOperator<numQubits>, m>;

	template<int numQubits>
	using MubSet = OperatorSet<numQubits, pow2(numQubits)>;

	template<int numQubits>
	using Mub = std::vector<MubSet<numQubits>>;


	template<int numQubits>
	auto expandStabilizer(const OperatorSet<numQubits, numQubits>& stabilizer) {
		MubSet<numQubits> expandedStabilizer{};
		for (int i = 0; i < pow2(numQubits); ++i) {
			for (int j = 0; j < numQubits; ++j) {
				if (i & (1 << j)) {
					expandedStabilizer[i].r += stabilizer[j].r;
					expandedStabilizer[i].s += stabilizer[j].s;
				}
			}
		}
		return expandedStabilizer;
	}

	template<int numQubits>
	constexpr void applyClifford(const Clifford<numQubits>& cliff, Mub<numQubits>& mub) {
		for (auto& set : mub) {
			for (auto& op : set) {
				op = op.applyClifford(cliff);
			}
		}
	}

	template<int numQubits>
	std::array<BinaryPauliOperator<numQubits>, numQubits> toEfficientStabilizer(const BinaryOperatorSet<numQubits, numQubits>& set) {
		std::array<BinaryPauliOperator<numQubits>, numQubits>eset;
		std::transform(set.begin(), set.end(), eset.begin(), [](const auto& op) {return BinaryPauliOperator<numQubits>{op}; });
		return eset;
	}

	template<int numQubits>
	MubSet<numQubits> toEfficientMubSet(const Q::MubSet<numQubits>& set) {
		MubSet<numQubits> eset;
		std::transform(set.begin(), set.end(), eset.begin(), [](const auto& op) {return BinaryPauliOperator<numQubits>{op}; });
		return eset;
	}

	template<int numQubits>
	Mub<numQubits> toEfficientMub(const Q::Mub<numQubits>& mub) {
		Mub<numQubits> emub(pow2(numQubits) + 1);
		std::transform(mub.begin(), mub.end(), emub.begin(), [](const auto& set) {return toEfficientMubSet(set); });
		return emub;
	}

	template<int numQubits, int m>
	Q::BinaryOperatorSet<numQubits, m> fromEfficientStabilizer(const OperatorSet<numQubits, m>& set) {
		Q::BinaryOperatorSet<numQubits, m> result;
		std::transform(set.begin(), set.end(), result.begin(), [](const auto& op) { return op.toBinaryPauliOperator(); });
		std::ranges::for_each(result, [](auto& op) {op.resetPhaseToTreatXZasY(); });
		return result;
	}

	template<int numQubits>
	Q::Mub<numQubits> fromEfficientMub(const Mub<numQubits>& emub) {
		Q::Mub<numQubits> mub(pow2(numQubits) + 1);
		for (size_t i = 0; i < pow2(numQubits) + 1; ++i) {
			std::transform(emub[i].begin(), emub[i].end(), mub[i].begin(), [](const auto& op) { return op.toBinaryPauliOperator(); });
			std::ranges::for_each(mub[i], [](auto& op) {op.resetPhaseToTreatXZasY(); });
		}
		return mub;
	}



	/// @brief Count operators that have identities at all positions specified by a given identity structure. 
	///        E.g. if identityStructure ist 011011, count operators in given set that are of the form II·II·
	///        with any of X,Y,Z at "·". 
	/// @param set Set of commuting Pauli operators
	template<int n>
	int countIdentityStructure(const MubSet<n>& set, const BinaryVector<n>& identityStructure) {
		return std::accumulate(set.begin(), set.end(), 0, [identityStructure](int count, auto op) {
			return count + (op.getIdentityString() == identityStructure); }
		);
	}

	/// @brief Get entry A_n of the sector length distribution associated with the given set. 
	/// @param set Set of commuting Pauli operators
	/// @return 
	template<int n>
	int getPauliWeightNCount(const MubSet<n>& set) {
		return countIdentityStructure(set, {});
	}

	/// @brief Returns a binary vector with 1 for every qubit that is entangled and 0 for every unentangled qubit. 
	template<int n>
	auto whichQubitsAreEntangled(const MubSet<n>& set) {
		BinaryVector<n> result;
		constexpr size_t entangledICount = std::max(1ULL, pow2(n - 2)) - 1;
		for (int i = 0; i < n; ++i) {
			int countI = std::accumulate(set.begin(), set.end(), 0ULL, [i](size_t count, const auto& op) { return count + op.getIdentityString().get(i); });
			result.set(i, countI == entangledICount);
		}
		return result;
	}



	/// @brief Get SLD characterization of a 2-qubit MUB
	///        Output index  |  SLD  |  graph
	///        --------------|-------|-------------
	///             [0]      | 1 0 3 | empty graph
	///             [1]      | 1 2 1 | one pair
	auto sldCharacterizeMub(const Mub<2>& mub) {
		std::array<int, 2> result{};
		for (const auto& set : mub) {
			int pauliWeight4Count = std::accumulate(set.begin(), set.end(), 0, [](int a, auto op) {
				//print("{},{}  ", op.toString(), op.getIdentityString());
				return a + (op.getIdentityString() == 0); }
			);
			//println("\n");
			switch (pauliWeight4Count) {
			case 3: ++result[0]; break;
			case 1: ++result[1]; break;
			default: assert(false);
			}
		}
		return result;
	}


	/// @brief Get SLD characterization of a 4-qubit MUB
	///        Output index  |  SLD      |  graph
	///        --------------|-----------|-------------
	///             [0]      | 1 4 6 4 1 | empty graph
	///             [1]      | 1 2 4 6 3 | one pair
	///             [2]      | 1 1 3 7 4 | one triple
	///             [3]      | 1 0 6 0 9 | two pairs or star
	///             [4]      | 1 0 2 8 5 | line
	auto sldCharacterizeMub(const Mub<4>& mub) {
		std::array<int, 5> result{};
		for (const auto& set : mub) {
			int pauliWeight4Count = getPauliWeightNCount(set);
			switch (pauliWeight4Count) {
			case 5: ++result[4]; break;
			case 9: ++result[3]; break;
			case 4: ++result[2]; break;
			case 3: ++result[1]; break;
			case 1: ++result[0]; break;
			default: assert(false);
			}
		}
		return result;
	}


	/// @brief Get SLD characterization of a 5-qubit MUB
	///        Output index  |  SLD            |  graph
	///        --------------|-----------------|-----------------------------
	///              [0]     | 1 5 10 10  5  1 | graph with no edges, 1+1+1+1 
	///              [1]     | 1 3  6 10  9  3 | one pair, 1+1+2               
	///              [2]     | 1 2  4 10 11  4 | one triple, 1+3              
	///              [3]     | 1 1  6  6  9  9 | two pairs 2+2 or 1+star          
	///              [4]     | 1 1  2 10 13  5 | 1+line graph                 
	///              [5]     | 1 0 10  0  5 16 | star                         
	///              [6]     | 1 0  6  4  9 12 | pair + triple                
	///              [7]     | 1 0  4  6 11 10 | T (pusteblume)               
	///              [8]     | 1 0  2  8 13  8 | line                         
	///              [9]     | 1 0  0 10 15  6 | cycle                        
	auto sldCharacterizeMub(const Mub<5>& mub) {
		std::array<int, 6> result{};
		for (const auto& set : mub) {
			int pauliWeight5Count = getPauliWeightNCount(set);
			switch (pauliWeight5Count) {
			case 6: ++result[9]; break;
			case 8: ++result[8]; break;
			case 10: ++result[7]; break;
			case 12: ++result[6]; break;
			case 16:++result[5]; break;
			case 5: ++result[4]; break;
			case 9: ++result[3]; break;
			case 4: ++result[2]; break;
			case 3: ++result[1]; break;
			case 1: ++result[0]; break;
			default: assert(false);
			}
		}
		return result;
	}

	template<int numQubits, class Generator>
	void randomClifford(Mub<numQubits>& mub, Generator& generator) {
		std::uniform_int_distribution<uint64_t> distribution{ 0, cliffordGroupSizeModuloPauli(numQubits) - 1 };
		auto symplectic = efficient::symplectic<numQubits>(distribution(generator));
		efficient::applyClifford(efficient::cliffordFrom2n2nSymplectic(symplectic), mub);
	}
}
