#pragma once
#include "binary_pauli.h"
#include "ht_circuits.h"
#include "special_math.h"
#include "formatting.h"
#include "pauli_operator_map.h"


namespace Q {

	template<int numQubits>
	using MubSet = BinaryOperatorSet<numQubits, pow2(numQubits) - 1>;

	template<int numQubits>
	using Mub = std::vector<MubSet<numQubits>>;




	auto printMUB = [](const auto& mub) {
		for (const auto& base : mub) {
			for (const auto& b : base) {
				std::cout << b.toString() << ' ';
			}
			std::cout << '\n';
		}
	};


	auto printMUBVertically = [](const auto& mub) {
		for (size_t i = 0; i < mub[0].size(); ++i) {
			for (size_t j = 0; j < mub.size(); ++j) {
				std::cout << mub[j][i].toString() << ' ';
			}
			std::cout << '\n';
		}
	};


	/// @brief Count operators with 0, 1, 2, ... , n occurrences of the Pauli operator op. 
	/// @return array with (numQubits+1) entries. The zeroth entry denotes for how many operators in the set
	///         op never occurs, the first entry how often op occurs once etc. 
	template<int numQubits>
	constexpr auto getCountStructure(const MubSet<numQubits>& mubSet, const BinaryPauliOperatorPrimitive& pauli) {
		std::array<int, numQubits + 1> counts{};
		for (const auto& oper : mubSet) {
			++counts[std::ranges::count(oper, pauli)];
		}
		return counts;
	}

	/// @brief Sector length distribution where the first number represents the number of 
	///        operators with Pauli weight 4, the second the number of operators with Pauli
	///        weight 3 and so on. 
	template<int numQubits>
	constexpr auto getSLD(const MubSet<numQubits>& mubSet) {
		auto sld = getCountStructure(mubSet, BinaryPauli::I);
		sld[numQubits] = 1; // a MUB set does not contain the I^n operator, so we add the entry for it here
		return sld;
	}

	/// @brief Count operators with 0, 1, 2, ... , n occurrences of the Pauli operator op. 
	/// @return array with (numQubits+1) entries. The zeroth entry denotes for how many operators in the set
	///         op never occurs, the first entry how often op occurs once etc. 
	template<int numQubits>
	constexpr auto getCountStructure(const Mub<numQubits>& mubSet, const BinaryPauliOperatorPrimitive& pauli) {
		std::array<std::array<int, numQubits + 1>, pow2(numQubits) + 1> counts;
		std::transform(mubSet.begin(), mubSet.end(), counts.begin(), [&pauli](const auto& mubSet) {return getCountStructure(mubSet, pauli); });
		return counts;
	}


	/// @brief Count the number of identity operators on the ith qubit
	template<int numQubits>
	int countPauliInMubSet(const MubSet<numQubits>& mubSet, int qubit, const BinaryPauliOperatorPrimitive& pauli) {
		return std::accumulate(mubSet.begin(), mubSet.end(), 0, [qubit, &pauli](int count, const auto& op) {return count + (op[qubit] == pauli); });
	}

	/// @brief For each set in the given mub count the number of identity operators on each qubit separately. 
	/// @return An array of arrays containing the counts. The inner array has numQubits entries, one for each qubit. 
	template<int numQubits>
	auto countPauliSetwise(const Mub<numQubits>& mub, const BinaryPauliOperatorPrimitive& pauli) {
		std::array<std::array<int, numQubits>, pow2(numQubits) + 1> counts{};
		for (size_t j = 0; j < mub.size(); ++j) {
			for (size_t i = 0; i < numQubits; ++i) {
				counts[j][i] = countPauliInMubSet(mub[j], i, pauli);
			}
		}
		return counts;
	}

	/// @brief Count the number of identity operators on the ith qubit
	template<int numQubits>
	int countIdentitiesInMubSet(const MubSet<numQubits>& mubSet, int qubit) {
		return countPauliInMubSet(mubSet, qubit, BinaryPauli::I);
	}


	/// @brief For each set in the given mub count the number of identity operators on each qubit separately. 
	/// @return An array of arrays containing the counts. The inner array has numQubits entries, one for each qubit. 
	template<int numQubits>
	auto countIdentitiesSetwise(const Mub<numQubits>& mub) {
		return countPauliSetwise(mub, BinaryPauli::I);
	}

	template<int n, int m>
	bool isQubitEntangled(const BinaryOperatorSet<n, m>& set, int qubit) {
		int qubitPauli = 0;
		for (int i = 0; i < m; ++i) {
			int pauli = (set[i].x(qubit).toInt() << 1) | set[i].z(qubit).toInt();
			if (pauli == 0) continue;
			if (qubitPauli == 0) {
				qubitPauli = pauli;
			}
			else {
				if (qubitPauli != pauli) return true;
			}
		}
		return false;
	}

	template<int n, int m>
	std::array<bool, n> areQubitsEntangled(const BinaryOperatorSet<n, m>& set) {
		std::array<bool, n> result;
		for (int i = 0; i < n; ++i) {
			result[i] = isQubitEntangled(set, i);
		}
		return result;
	}


	template<int numQubits>
	void swapQubitsInMub(Mub<numQubits>& mub, int qubit1, int qubit2) {
		std::ranges::for_each(mub, [&](auto& set) { std::ranges::for_each(set, [&](auto& op) { std::swap(op[qubit1], op[qubit2]); }); });
	}

	template<int numQubits>
	bool areMubwiseCommuting(const Mub<numQubits>& mub) {
		for (const auto& mubSet : mub) {
			for (int i = 0; i < mubSet.size(); ++i) {
				for (int j = 0; j < mubSet.size(); ++j) {
					if (commutator(mubSet[i], mubSet[j]) == 1)
						return false;
				}
			}
		}
		return true;
	}

	template<int n>
	bool isMub(const Mub<n>& mub) {
		PauliOperatorMap<int, n> map;
		map[PauliIndex<n>{ 0 }] = 1;
		for (const auto& mubSet : mub) {
			for (const auto& op : mubSet) {
				map[PauliIndex<n>{ op }] = 1;
			}
		}
		for (auto z : map) {
			if (z == 0) return false;
		}
		//return true;
		return areMubwiseCommuting(mub);
	}



	/// @brief Given one MUB set with 2^n-1 n-qubit operators find the "canonical" generating set 
	///        with respect to a diagonalization circuit. Here, canonical means that the diagonalized 
	///        operators have the form ZIII.., IZII..,IIZI,.. (possibly with a minus sign). The operators
	///        in the generating set are by construction independant. 
	///        
	/// @param mubSet                  A mutually unbiased basis. 
	/// @param diagonalizationCircuit  Circuit to use for the diagonalization
	/// @return Canonical generating set in the order that their diagonalized operators are ZIII.., IZII..,IIZI,..
	template<int numQubits>
	auto findCanonicalGeneratingSet(const MubSet<numQubits>& mubSet, const HTCircuit<numQubits>& diagonalizationCircuit) {
		BinaryOperatorSet<numQubits, numQubits> generatingSet;
		for (const auto& op : mubSet) {
			const auto resultOp = diagonalizationCircuit.transformPauli(op);
			if (std::ranges::count(resultOp, BinaryPauli::Z) == 1) {
				if (auto it = std::ranges::find(resultOp, BinaryPauli::Z); it != resultOp.end()) {
					generatingSet[std::distance(resultOp.begin(), it)] = op;
				}
				else {
					assert(false);
				}
			}
		}
		return generatingSet;
	}



	/// @brief Construct a MUB set with 2^n-1 operators from a set of n operators which form a 
	///        "canonical" generating set with respect to a diagonalization circuit. The input is
	///        is expected to be diagonalized by such a circuit to ZIII.., IZII..,IIZI,.. (possibly 
	///        with a minus sign). 
	/// 
	///        This is basically the inverse operation to findCanonicalGeneratingSet(). 
	///        
	/// @param generatingSet  A canonical generating set in the order that their diagonalized operators are ZIII.., IZII..,IIZI,..
	/// @return A complete MUB set with 2^n-1 operators
	template<int numQubits>
	auto constructMubSetFromCanonicalGeneratingSet(const BinaryOperatorSet<numQubits, numQubits>& generatingSet) {
		MubSet<numQubits> mubSet;
		for (size_t i = 1; i < pow2(numQubits); ++i) {
			for (size_t j = 0; j < numQubits; ++j) {
				if (i & (1ULL << j))
					mubSet[i - 1] *= generatingSet[j];
			}
			mubSet[i - 1].resetPhaseToTreatXZasY();
		}
		return mubSet;
	}


	/// @brief Write a string of a form like
	///               XZI YXX IZY : n=3:403:0-1,1-2,
	///        containing the canonical generating set and a diagonalization circuit for the 
	///        given mutually unbiased basis. 
	/// 
	/// @param mubSet                  A mutually unbiased basis. 
	/// @param diagonalizationCircuit  Circuit to use for the diagonalization
	/// @return a string
	template<int numQubits>
	auto writeMUBDiagonalizationSpecification(const MubSet<numQubits>& mubSet, const HTCircuit<numQubits>& diagonalizationCircuit) {
		//const auto generatingSet = findCanonicalGeneratingSet(mubSet, diagonalizationCircuit);
		const auto generatingSet = diagonalizationCircuit.toQuantumCircuit().inverse().getStabilizer();
		std::string result;
		for (const auto& op : generatingSet) {
			result += std::format("{} ", op);
		}
		result += ": ";
		result += diagonalizationCircuit.serialize();
		return result;
	}



}
