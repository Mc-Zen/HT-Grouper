
#pragma once
#include <cstdint>
#include <string_view>
#include "binary_pauli.h"
#include <format>


namespace Q {

	struct Pauli {
	public:
		uint64_t r{};
		uint64_t s{};
		int n{ 1 };
		BinaryPhase phase;

		constexpr Pauli() = default;

		/// @brief Creates an identity Pauli operator of length n
		/// @param n Number of qubits
		explicit constexpr Pauli(int n) : n(n) {};

		explicit(false) constexpr Pauli(std::string_view pauliString) {
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

		static constexpr Pauli SingleX(int n, int qubit) {
			Pauli pauli{ n };
			pauli.setX(qubit, 1);
			return pauli;
		}

		static constexpr Pauli SingleZ(int n, int qubit) {
			Pauli pauli{ n };
			pauli.setZ(qubit, 1);
			return pauli;
		}


		constexpr auto x(int qubit) const { return (r >> qubit) & 1ULL; }
		constexpr auto z(int qubit) const { return (s >> qubit) & 1ULL; }
		constexpr void setX(int qubit, bool value) { r ^= (-value ^ r) & (1ULL << qubit); }
		constexpr void setZ(int qubit, bool value) { s ^= (-value ^ r) & (1ULL << qubit); }


		/// @brief Get phase of the operator like when XZ is represented as -iY
		constexpr BinaryPhase getPhase() const { return phase - getYPhase(); }

		/// @brief Get phase of the operator like when Y is represented as iXZ
		constexpr BinaryPhase getXZPhase() const { return phase; }

		constexpr void increasePhase(int phaseInc) { phase += phaseInc; }
		constexpr void decreasePhase(int phaseDec) { phase -= phaseDec; }

		int pauliWeight() const { return std::popcount(r | s); }
		int identityCount() const { return n - pauliWeight(); }
		auto getIdentityString() const { return ~(r | s); }

		std::string toString() const {
			constexpr std::array<char, 4> c{ 'I','X','Z','Y' };
			std::string str;
			str.reserve(n);
			for (int i = 0; i < n; ++i) str += c[x(i) + z(i) * 2];
			return str;
		}

		constexpr friend bool operator==(const Pauli& a, const Pauli& b) = default;

		/// @brief Commutator of two Pauli operators. The result is in binary form, 0 if 
		///        p1 and p2 commute, 1 if they anticommute. 
		constexpr friend int commutator(const Pauli& p1, const Pauli& p2) {
			return std::popcount((p1.r & p2.s) ^ (p2.r & p1.s)) & 1;
		}

	private:
		constexpr void fromStringOperator(const std::string_view& str) {
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

		// Get the phase that is accumulated by representing Y as iXZ
		constexpr BinaryPhase getYPhase() const { return { std::popcount(r & s) }; }

	};

}


template<class CharT>
struct std::formatter<Q::Pauli, CharT> : std::formatter<std::string_view, CharT> {
	template<class FormatContext>
	auto format(const Q::Pauli& op, FormatContext& fc) const {
		const auto phase = op.getPhase();
		if (phase != Q::BinaryPhase{ 0 })
			std::format_to(fc.out(), "{}", phase.toString());
		return std::format_to(fc.out(), "{}", op.toString());
	}
};