
#pragma once
#include "matrix.h"
#include "binary.h"
#include "binary_phase.h"
#include <string>
#include <algorithm>

namespace Q {

	template<int m, int n>
	using BinaryMatrix = Math::Matrix<Binary, m, n>;

	template<int n>
	using BinaryVector = Math::Vector<Binary, n>;


	// single qubit pauli operator in the form X^r Z^s with myBinaryPauli[0] = r, myBinaryPauli[1] = s
	using BinaryPauliOperatorPrimitive = BinaryVector<2>;
	// single qubit bit clifford gate in binary form
	using BinaryCliffordGate = BinaryMatrix<2, 2>;


	namespace BinaryPauli {
		inline constexpr auto I = BinaryPauliOperatorPrimitive{ {0, 0} };
		inline constexpr auto X = BinaryPauliOperatorPrimitive{ {1, 0} };
		inline constexpr auto Y = BinaryPauliOperatorPrimitive{ {1, 1} };
		inline constexpr auto Z = BinaryPauliOperatorPrimitive{ {0, 1} };
	}

	namespace BinaryCliffordGates {
		inline constexpr auto I = BinaryCliffordGate{ 1,0,0,1 };
		inline constexpr auto H = BinaryCliffordGate{ 0,1,1,0 };
		inline constexpr auto S = BinaryCliffordGate{ 1,0,1,1 };
		inline constexpr auto HS = BinaryCliffordGate{ 1,1,1,0 };
		inline constexpr auto SH = BinaryCliffordGate{ 0,1,1,1 };
		inline constexpr auto HSH = BinaryCliffordGate{ 1,1,0,1 };
	}

	constexpr char toChar(const BinaryPauliOperatorPrimitive& op) {
		constexpr std::array<char, 4> c{ 'I','X','Z','Y' };
		return c[op[0].toInt() + op[1].toInt() * 2];
	}

	constexpr BinaryPauliOperatorPrimitive binaryPauliOperatorPrimitiveFromChar(char c) {
		switch (c) {
		case 'I': return BinaryPauli::I;
		case 'X': return BinaryPauli::X;
		case 'Y': return BinaryPauli::Y;
		case 'Z': return BinaryPauli::Z;
		default: return BinaryPauli::I;
		}
	}



	/// @brief Binary n-qubit pauli operator with phase in the form i^q with q=0,1,2,3
	template<int n>
	class BinaryPauliOperator {
	public:
		std::array<BinaryPauliOperatorPrimitive, n> ops;
		BinaryPhase phase;

		constexpr BinaryPauliOperator() = default;

		// Accepted inputs: XIX, XYZ, +XXI, -IXY, -iXXX, iZZZ. 
		explicit(false) constexpr BinaryPauliOperator(const std::string_view& sv);


		/// @brief Generate a Pauli operator of the form X^r (e.g. XIIIXX)
		/// @param r Bit string with LSB for first qubit
		static constexpr BinaryPauliOperator FromXString(uint64_t r);

		/// @brief Generate a Pauli operator of the form Z^s (e.g. ZIIZIZZ)
		/// @param s Bit string with LSB for first qubit
		static constexpr BinaryPauliOperator FromZString(uint64_t s);

		/// @brief Create a Binary pauli operator that has Z at given index
		static constexpr BinaryPauliOperator SingleZ(int index) { BinaryPauliOperator op; op[index] = BinaryPauli::Z; return op; }

		/// @brief Create a Binary pauli operator that has X at given index
		static constexpr BinaryPauliOperator SingleX(int index) { BinaryPauliOperator op; op[index] = BinaryPauli::X; return op; }

		/// @brief Get phase of the operator like when XZ is represented as -iY
		constexpr BinaryPhase getPhase() const { return phase - getYPhase(); }

		/// @brief Get phase of the operator like when Y is represented as iXZ
		constexpr BinaryPhase getXZPhase() const { return phase; }

		constexpr void increasePhase(int phaseInc) { phase += phaseInc; }
		constexpr void decreasePhase(int phaseDec) { phase -= phaseDec; }


		constexpr BinaryPauliOperatorPrimitive& operator[](size_t i) { return ops[i]; }
		constexpr const BinaryPauliOperatorPrimitive& operator[](size_t i) const { return ops[i]; }

		constexpr Binary& x(size_t i) { return ops[i][0]; }
		constexpr const Binary& x(size_t i) const { return ops[i][0]; }
		constexpr Binary& z(size_t i) { return ops[i][1]; }
		constexpr const Binary& z(size_t i) const { return ops[i][1]; }

		constexpr auto begin() { return ops.begin(); }
		constexpr auto end() { return ops.end(); }
		constexpr auto begin() const { return ops.begin(); }
		constexpr auto end() const { return ops.end(); }
		constexpr auto cbegin() const { return begin(); }
		constexpr auto cend() const { return end(); }

		std::string toString(bool printPhase = false) const;

		/// @brief BinaryPauliOperator describes an operator i^qX^rZ^s. Get r with first qubit at LSB
		constexpr uint64_t getXString() const;

		/// @brief BinaryPauliOperator describes an operator i^qX^rZ^s. Get s with first qubit at LSB
		constexpr uint64_t getZString() const;

		constexpr int identityCount() const { return std::ranges::count(ops, BinaryPauli::I); }
		constexpr int pauliWeight() const { return n - identityCount(); }

		// Reset the internal phase to match the number of Y operators (used to represent Y operators in XZ form)
		constexpr void resetPhaseToTreatXZasY() {
			phase = BinaryPhase{ 0 };
			for (const auto& op : ops) phase += op[0] & op[1];
		}

		// Apply other binary pauli operator
		constexpr BinaryPauliOperator& operator*=(const BinaryPauliOperator& other);

		constexpr friend bool operator==(const BinaryPauliOperator& a, const BinaryPauliOperator& b) = default;
		constexpr friend bool operator==(const BinaryPauliOperator& a, const std::string_view& b) { return a == BinaryPauliOperator{ b }; }


	private:
		constexpr void fromStringOperator(const std::string_view& str) {
			std::transform(str.begin(), str.end(), ops.begin(), &binaryPauliOperatorPrimitiveFromChar);
		}
		// Get the phase that is accumulated by representing Y as iXZ
		constexpr BinaryPhase getYPhase() const {
			BinaryPhase yPhase{ 0 };
			for (const auto& op : ops) yPhase += op[0] & op[1];
			return yPhase;
		}

	};


	namespace Clifford {
		template<int n>
		void x(BinaryPauliOperator<n>& pauli, int qubit) {
			pauli.increasePhase(2 * pauli.z(qubit).toInt());
		}

		template<int n>
		void y(BinaryPauliOperator<n>& pauli, int qubit) {
			pauli.increasePhase(2 * (pauli.x(qubit) + pauli.z(qubit)).toInt());
		}

		template<int n>
		void z(BinaryPauliOperator<n>& pauli, int qubit) {
			pauli.increasePhase(2 * pauli.x(qubit).toInt());
		}

		template<int n>
		void h(BinaryPauliOperator<n>& pauli, int qubit) {
			std::swap(pauli.x(qubit), pauli.z(qubit));
			pauli.increasePhase(2 * (pauli.x(qubit) * pauli.z(qubit)).toInt());
		}

		template<int n>
		void s(BinaryPauliOperator<n>& pauli, int qubit) {
			pauli.z(qubit) += pauli.x(qubit);
			pauli.increasePhase(pauli.x(qubit).toInt());
		}

		template<int n>
		void sdg(BinaryPauliOperator<n>& pauli, int qubit) {
			pauli.z(qubit) += pauli.x(qubit);
			pauli.decreasePhase(pauli.x(qubit).toInt());
		}

		template<int n>
		void hs(BinaryPauliOperator<n>& pauli, int qubit) {
			s(pauli, qubit);
			h(pauli, qubit);
		}

		template<int n>
		void sh(BinaryPauliOperator<n>& pauli, int qubit) {
			h(pauli, qubit);
			s(pauli, qubit);
		}

		template<int n>
		void hsh(BinaryPauliOperator<n>& pauli, int qubit) {
			h(pauli, qubit);
			s(pauli, qubit);
			h(pauli, qubit);
		}

		template<int n>
		void cx(BinaryPauliOperator<n>& pauli, int control, int target) {
			pauli.x(target) += pauli.x(control);
			pauli.z(control) += pauli.z(target);
		}

		template<int n>
		void cz(BinaryPauliOperator<n>& pauli, int qubit1, int qubit2) {
			pauli.z(qubit2) += pauli.x(qubit1);
			pauli.z(qubit1) += pauli.x(qubit2);
			pauli.increasePhase(2 * (pauli.x(qubit1) * pauli.x(qubit2)).toInt()); // if both operators have X component: phase flip
		}

		template<int n>
		void swap(BinaryPauliOperator<n>& pauli, int qubit1, int qubit2) {
			std::swap(pauli[qubit1], pauli[qubit2]);
		}
	}


	template<int n>
	constexpr BinaryPauliOperator<n>::BinaryPauliOperator(const std::string_view& sv) {
		assert(sv.size() >= n);
		if (sv.size() > n) {
			if (sv.starts_with('i'))
				phase += 1;
			else if (sv.starts_with("-i")) phase += 3;
			else if (sv.starts_with('-')) phase += 2;
			fromStringOperator(sv.substr(sv.size() - n, n));
		}
		else {
			fromStringOperator(sv);
		}
		phase += getYPhase();
	}

	template<int n>
	static constexpr BinaryPauliOperator<n> BinaryPauliOperator<n>::FromXString(uint64_t r) {
		BinaryPauliOperator op;
		for (int j = 0; j < n; ++j) {
			if (r & (1ULL << j)) {
				op[j] = BinaryPauli::X;
			}
		}
		return op;
	}

	template<int n>
	static constexpr BinaryPauliOperator<n> BinaryPauliOperator<n>::FromZString(uint64_t s) {
		BinaryPauliOperator op;
		for (int j = 0; j < n; ++j) {
			if (s & (1ULL << j)) {
				op[j] = BinaryPauli::Z;
			}
		}
		return op;
	}

	template<int n>
	constexpr BinaryPauliOperator<n>& BinaryPauliOperator<n>::operator*=(const BinaryPauliOperator<n>& other) {
		for (size_t i = 0; i < n; i++) {
			ops[i] += other.ops[i];
		}
		phase += other.phase;
		return *this;
	}

	template<int n>
	std::string BinaryPauliOperator<n>::toString(bool printPhase) const {
		std::string s;
		s.reserve(n + 1);
		if (printPhase)
			s += phase.toString();
		for (const auto& op : ops) s += toChar(op);
		return s;
	}

	template<int n>
	constexpr uint64_t BinaryPauliOperator<n>::getXString() const {
		uint64_t xString{};
		for (int i = 0; i < n; ++i) {
			xString <<= 1;
			xString |= ops[n - i - 1][0];
		}
		return xString;
	}

	template<int n>
	constexpr uint64_t BinaryPauliOperator<n>::getZString() const {
		uint64_t zString{};
		for (int i = 0; i < n; ++i) {
			zString <<= 1;
			zString |= ops[n - i - 1][1];
		}
		return zString;
	}




	template<int numQubits, int numOperators>
	using BinaryOperatorSet = std::array<BinaryPauliOperator<numQubits>, numOperators>;





	namespace MUBTransforms {

		constexpr void localXZSwap(BinaryPauliOperatorPrimitive& op) {
			std::swap(op[0], op[1]);
		}

		constexpr void localXYSwap(BinaryPauliOperatorPrimitive& op) {
			if (op[0]) op[1] += 1;
		}

		constexpr void localYZSwap(BinaryPauliOperatorPrimitive& op) {
			if (op[1]) op[0] += 1;
		}

		constexpr void localPermutationXYZ(BinaryPauliOperatorPrimitive& op) {
			int a = (op[1].toInt() << 1) + op[0].toInt(); // 0:I, 1:X, 2:Z, 3:Y
			if (a == 0) return;
			a = (a + 1) % 3 + 1;
			op[0] = a & 0b01;
			op[1] = a & 0b10;
		}

		template<int n>
		constexpr void localXZSwap(BinaryPauliOperator<n>& op, int qubit) {
			localXZSwap(op.ops[qubit]);
		}

		template<int n>
		constexpr void localXYSwap(BinaryPauliOperator<n>& op, int qubit) {
			localXYSwap(op.ops[qubit]);
			op.resetPhaseToTreatXZasY();
		}

		template<int n>
		constexpr void localYZSwap(BinaryPauliOperator<n>& op, int qubit) {
			localYZSwap(op.ops[qubit]);
			op.resetPhaseToTreatXZasY();
		}

		template<int n>
		constexpr void localPermutationXYZ(BinaryPauliOperator<n>& op, int qubit) {
			localPermutationXYZ(op.ops[qubit]);
			op.resetPhaseToTreatXZasY();
		}


		template<int n, int m>
		constexpr void localXZSwap(std::array<BinaryPauliOperator<n>, m>& ops, int qubit) {
			for (auto& op : ops) localXZSwap(op, qubit);
		}

		template<int n, int m>
		constexpr void localXYSwap(std::array<BinaryPauliOperator<n>, m>& ops, int qubit) {
			for (auto& op : ops) localXYSwap(op, qubit);
		}

		template<int n, int m>
		constexpr void localYZSwap(std::array<BinaryPauliOperator<n>, m>& ops, int qubit) {
			for (auto& op : ops) localYZSwap(op, qubit);
		}

		template<int n, int m>
		constexpr void localPermutationXYZ(std::array<BinaryPauliOperator<n>, m>& ops, int qubit) {
			for (auto& op : ops) localPermutationXYZ(op, qubit);
		}



		constexpr void applyCZ(BinaryPauliOperatorPrimitive& op1, BinaryPauliOperatorPrimitive& op2) {
			op1[1] += op2[0];
			op2[1] += op1[0];
		}

		template<int n>
		constexpr void applyCZ(BinaryPauliOperator<n>& op, int qubit1, int qubit2) {
			applyCZ(op[qubit1], op[qubit2]);
			op.phase += 2 * (op[qubit1][0] & op[qubit2][0]); // if both operators have X component: phase flip
		}

		template<int n, int m>
		constexpr void applyCZ(std::array<BinaryPauliOperator<n>, m>& ops, int qubit1, int qubit2) {
			for (auto& op : ops) applyCZ(op, qubit1, qubit2);
		}

		constexpr void applyCX(BinaryPauliOperatorPrimitive& control, BinaryPauliOperatorPrimitive& target) {
			target[0] += control[0];
			control[1] += target[1];
		}

		template<int n>
		constexpr void applyCX(BinaryPauliOperator<n>& op, int control, int target) {
			//if (op[target] == BinaryPauli::Y) op.phase += 2;
	/*		if ((op[control] == BinaryPauli::Y && op[target] == BinaryPauli::Y) || (op[target] == BinaryPauli::Z && op[control] == BinaryPauli::X)) {
				op.phase += 2;
			}*/
			//op.phase += 2 * (op[target][1] & op[control][0]); // if both operators have X component: phase flip
			applyCX(op[control], op[target]);
			//if (op[target] == BinaryPauli::Y) op.phase += 2;


		}

		template<int n, int m>
		constexpr void applyCX(std::array<BinaryPauliOperator<n>, m>& ops, int control, int target) {
			for (auto& op : ops) applyCX(op, control, target);
		}



	}







	constexpr std::string toString(const BinaryCliffordGate& g) {
		if (g == BinaryCliffordGates::I) return "I";
		if (g == BinaryCliffordGates::H) return "H";
		if (g == BinaryCliffordGates::S) return "S";
		if (g == BinaryCliffordGates::SH) return "SH";
		if (g == BinaryCliffordGates::HSH) return "HSH";
		if (g == BinaryCliffordGates::HS) return "HS";
		return "";
	}




	constexpr Binary commutator(const BinaryPauliOperatorPrimitive& b1, const BinaryPauliOperatorPrimitive& b2) {
		return b1[0] * b2[1] + b1[1] * b2[0];
	}

	template<int n>
	constexpr Binary commutator(const BinaryPauliOperator<n>& b1, const BinaryPauliOperator<n>& b2) {
		Binary c;
		for (size_t i = 0; i < n; ++i) c += commutator(b1[i], b2[i]);
		return c;
	}


	/// @brief generate a MUB set from a string of operators, e.g. parseMubSet<3,4>("XYX XXY ZXY ZZZ")
	/// @param string input sequence
	/// @return MUB
	template<int n, int m>
	constexpr std::array<BinaryPauliOperator<n>, m> parseMubSet(const std::string_view& string) {
		std::array<BinaryPauliOperator<n>, m> result{};
		assert(string.size() == m * n + (m - 1));

		for (size_t i = 0; i < m; i++) {
			result[i] = { string.substr((n + 1) * i,n) };
		}
		return result;
	}

}


