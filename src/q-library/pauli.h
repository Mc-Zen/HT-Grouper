
#pragma once
#include <cstdint>
#include <string_view>
#include "formatting.h"
#include "binary_phase.h"


namespace Q {

	/// @brief Representation of a Pauli operator. 
	struct Pauli {
	public:
		using Bitstring = uint64_t;

		constexpr Pauli() = default;

		/// @brief Creates an identity Pauli operator of length n
		/// @param n Number of qubits
		explicit constexpr Pauli(int n) : n(n) {};

		/// @brief Create a Pauli operator from a string, e.g. XIIXZ, -XYYYX, -iZZ, iXIX
		explicit(false) constexpr Pauli(std::string_view pauliString);

		/// @brief Create A Pauli operator with just one X at the specified position, e.g. IIXIII
		static constexpr Pauli SingleX(int n, int qubit);
		/// @brief Create A Pauli operator with just one Z at the specified position, e.g. IIZIII
		static constexpr Pauli SingleZ(int n, int qubit);

		static constexpr Pauli Identity(int n);



		constexpr int numQubits() const { return n; };

		/// @brief Check if operator at given qubit has an x component
		constexpr uint64_t x(int qubit) const;
		/// @brief Check if operator at given qubit has a z component
		constexpr uint64_t z(int qubit) const;

		/// @brief Set x component at given qubit position. The behaviour is unspecified if value is not 0 or 1
		constexpr void setX(int qubit, int value);
		/// @brief Set z component at given qubit position. The behaviour is unspecified if value is not 0 or 1
		constexpr void setZ(int qubit, int value);


		/// @brief Get phase of the operator when XZ is represented as -iY
		constexpr BinaryPhase getPhase() const { return phase - getYPhase(); }

		/// @brief Get phase of the operator when Y is represented as iXZ
		constexpr BinaryPhase getXZPhase() const { return phase; }

		constexpr void increasePhase(int phaseInc) { phase += phaseInc; }
		constexpr void decreasePhase(int phaseDec) { phase -= phaseDec; }

		/// @brief Get the operators Pauli weight (i.e. the number of non-identity single-qubit Paulis)
		constexpr int pauliWeight() const;
		/// @brief Get the number of identity single-qubit Paulis in the operator
		constexpr int identityCount() const;

		/// @brief Get the X components of the Pauli as binary string, e.g. XYZI -> 1100
		constexpr Bitstring getXString() const;
		/// @brief Get the Z components of the Pauli as binary string, e.g. XYZI -> 0110
		constexpr Bitstring getZString() const;

		/// @brief Get a binary string with 1 for each identity, e.g. XYZI -> 0001
		constexpr Bitstring getIdentityString() const;

		constexpr std::string toString() const;

		constexpr friend bool operator==(const Pauli& a, const Pauli& b) = default;

		/// @brief Commutator of two Pauli operators. The result is in binary form, 0 if 
		///        p1 and p2 commute, 1 if they anticommute. 
		constexpr friend int commutator(const Pauli& p1, const Pauli& p2);

		/// @brief Check if p1 and p2 commute on each qubit. 
		constexpr friend bool commutesQubitWise(const Pauli& p1, const Pauli& p2);

		/// @brief Check if p1 and p2 commute locally on subset A, i.e. whether p1' and p2' commute
		///        where 
		///            p'[i] = | p[i]  if i in A
		///                    | I     else.
		///        The argument support encodes A with a bit at location j set to 1 if j in A. 
		constexpr friend bool commutesLocally(const Pauli& p1, const Pauli& p2, uint64_t support);


	private:
		constexpr void fromStringOperator(const std::string_view& str);

		// Get the phase that is accumulated by representing Y as iXZ
		constexpr BinaryPhase getYPhase() const { return { std::popcount(r & s) }; }

		Bitstring r{};
		Bitstring s{};
		int n{ 1 };
		BinaryPhase phase;
	};





	constexpr Pauli::Pauli(std::string_view pauliString) {
		if (pauliString.starts_with('i')) {
			phase += 1;
			fromStringOperator(pauliString.substr(1));
		}
		else if (pauliString.starts_with("-i")) {
			phase += 3;
			fromStringOperator(pauliString.substr(2));
		}
		else if (pauliString.starts_with('-')) {
			phase += 2;
			fromStringOperator(pauliString.substr(1));
		}
		else {
			fromStringOperator(pauliString);
		}
		phase += getYPhase();
	}

	constexpr Pauli Pauli::SingleX(int n, int qubit) {
		Pauli pauli{ n };
		pauli.setX(qubit, 1);
		return pauli;
	}

	constexpr Pauli Pauli::SingleZ(int n, int qubit) {
		Pauli pauli{ n };
		pauli.setZ(qubit, 1);
		return pauli;
	}

	constexpr Pauli Pauli::Identity(int n) {
		return Pauli{ n };
	}


	constexpr uint64_t Pauli::x(int qubit) const { return (r >> qubit) & 1ULL; }

	constexpr uint64_t Pauli::z(int qubit) const { return (s >> qubit) & 1ULL; }

	constexpr void Pauli::setX(int qubit, int value) { r ^= (-value ^ r) & (1ULL << qubit); }

	constexpr void Pauli::setZ(int qubit, int value) { s ^= (-value ^ s) & (1ULL << qubit); }

	/// @brief Get phase of the operator like when XZ is represented as -iY

	constexpr int Pauli::pauliWeight() const { return std::popcount(r | s); }

	constexpr int Pauli::identityCount() const { return n - pauliWeight(); }

	constexpr Pauli::Bitstring Pauli::getIdentityString() const { return ~(r | s); }

	constexpr Pauli::Bitstring Pauli::getXString() const { return r; }

	constexpr Pauli::Bitstring Pauli::getZString() const { return s; }

	constexpr std::string Pauli::toString() const {
		constexpr std::array<char, 4> c{ 'I','X','Z','Y' };
		std::string str;
		str.reserve(n);
		for (int i = 0; i < n; ++i) str += c[x(i) + z(i) * 2];
		return str;
	}


	constexpr void Pauli::fromStringOperator(const std::string_view& str) {
		n = static_cast<int>(str.length());
		int i{};
		for (char c : str) {
			switch (c) {
			case 'I': break;
			case 'X': r |= (1ULL << i); break;
			case 'Y': r |= (1ULL << i); s |= (1ULL << i); break;
			case 'Z': s |= (1ULL << i); break;
			default: break;
			}
			++i;
		}
	}

	/// @brief Commutator of two Pauli operators. The result is in binary form, 0 if 
	///        p1 and p2 commute, 1 if they anticommute. 
	constexpr int commutator(const Pauli& p1, const Pauli& p2) { return std::popcount((p1.r & p2.s) ^ (p2.r & p1.s)) & 1; }


	constexpr bool commutesQubitWise(const Pauli& p1, const Pauli& p2) { return ~(p1.getIdentityString() | p2.getIdentityString() | (~(p1.getXString() ^ p2.getXString()) & ~(p1.getZString() ^ p2.getZString()))) == 0; }


	constexpr bool commutesLocally(const Pauli& p1, const Pauli& p2, uint64_t support) {
		return (std::popcount((p1.r & p2.s & support) ^ (p2.r & p1.s & support)) & 1) == 0;
	}

	///// @brief See @commutesLocally(const Pauli& p1, const Pauli& p2, int64_t support), 
	/////        but here the A is directly stored as indices in a container. 
	//template<class ForwardIterable> requires requires (ForwardIterable container) { {std::begin(container) } -> std::convertible_to<typename ForwardIterable::iterator>; }
	//constexpr bool commutesLocally(const Pauli& p1, const Pauli& p2, const ForwardIterable& A);

	//template<class ForwardIterable>
	//constexpr bool commutesLocally(const Pauli& p1, const Pauli& p2, const ForwardIterable& A) {
	//	int64_t support{};
	//	for (auto index : A) {
	//		support |= (1ULL << index);
	//	}
	//	return commutesLocally(p1, p2, support);
	//}






	namespace Clifford {
		constexpr void x(Pauli& pauli, int qubit) {
			pauli.increasePhase(2 * pauli.z(qubit));
		}

		constexpr void y(Pauli& pauli, int qubit) {
			pauli.increasePhase(2 * (pauli.x(qubit) + pauli.z(qubit)));
		}

		constexpr void z(Pauli& pauli, int qubit) {
			pauli.increasePhase(2 * pauli.x(qubit));
		}

		constexpr void h(Pauli& pauli, int qubit) {
			auto x = pauli.x(qubit);
			auto z = pauli.z(qubit);
			pauli.setX(qubit, z);
			pauli.setZ(qubit, x);
			pauli.increasePhase(2 * (pauli.x(qubit) * pauli.z(qubit)));
		}

		constexpr void s(Pauli& pauli, int qubit) {
			pauli.setZ(qubit, pauli.z(qubit) ^ pauli.x(qubit));
			pauli.increasePhase(pauli.x(qubit));
		}

		constexpr void sdg(Pauli& pauli, int qubit) {
			pauli.setZ(qubit, pauli.z(qubit) ^ pauli.x(qubit));
			pauli.decreasePhase(pauli.x(qubit));
		}

		constexpr void hs(Pauli& pauli, int qubit) {
			s(pauli, qubit);
			h(pauli, qubit);
		}

		constexpr void sh(Pauli& pauli, int qubit) {
			h(pauli, qubit);
			s(pauli, qubit);
		}

		constexpr void hsh(Pauli& pauli, int qubit) {
			h(pauli, qubit);
			s(pauli, qubit);
			h(pauli, qubit);
		}

		constexpr void cx(Pauli& pauli, int control, int target) {
			pauli.setX(target, pauli.x(target) ^ pauli.x(control));
			pauli.setZ(control, pauli.z(control) ^ pauli.z(target));
		}

		constexpr void cz(Pauli& pauli, int qubit1, int qubit2) {
			pauli.setZ(qubit2, pauli.z(qubit2) ^ pauli.x(qubit1));
			pauli.setZ(qubit1, pauli.z(qubit1) ^ pauli.x(qubit2));
			pauli.increasePhase(2 * (pauli.x(qubit1) * pauli.x(qubit2))); // if both operators have X component: phase flip
		}

		constexpr void swap(Pauli& pauli, int qubit1, int qubit2) {
			auto x1 = pauli.x(qubit1);
			auto z1 = pauli.z(qubit1);
			auto x2 = pauli.x(qubit2);
			auto z2 = pauli.z(qubit2);
			pauli.setX(qubit1, x2);
			pauli.setZ(qubit1, z2);
			pauli.setX(qubit2, x1);
			pauli.setZ(qubit2, z1);
		}
	}
}


template<class CharT>
struct std::formatter<Q::Pauli, CharT> : std::formatter<std::string_view, CharT> {
	template<class FormatContext>
	auto format(const Q::Pauli& op, FormatContext& fc) const {
		if (const auto phase = op.getPhase(); phase != Q::BinaryPhase{ 0 }) {
			std::format_to(fc.out(), "{}", phase.toString());
		}
		return std::format_to(fc.out(), "{}", op.toString());
	}
};