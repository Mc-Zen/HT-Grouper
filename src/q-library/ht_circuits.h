#pragma once

#include "binary_pauli.h"
#include "special_math.h"
#include "quantum_circuit.h"
#include "graph.h"
#include <bit>
#include <vector>
#include <ranges>
#include <iostream>

namespace Q {


	template<int numQubits = 2>
	class HTCircuit {
	public:
		using AdjencencyMatrix = BinaryMatrix<numQubits, numQubits>;

		Graph<numQubits> graph;
		std::array<BinaryCliffordGate, numQubits> singleQubitLayer;

		auto transformPauli(const BinaryPauliOperator<numQubits>& input) const {
			BinaryPauliOperator<numQubits> op{ input };
			transformThroughSingleQubitLayer(op);
			transformThroughCZ(op);
			transformThroughHadamardLayer(op);
			return op;
		}


		auto toQuantumCircuit() const {
			QuantumCircuit<numQubits> qc;
			size_t k = 0;
			for (const auto& gate : singleQubitLayer) {
				if (gate == BinaryCliffordGates::H) qc.h(k);
				else if (gate == BinaryCliffordGates::S) qc.s(k);
				else if (gate == BinaryCliffordGates::SH) { qc.h(k); qc.s(k); }
				else if (gate == BinaryCliffordGates::HSH) { qc.h(k); qc.s(k); qc.h(k); }
				else if (gate == BinaryCliffordGates::HS) { qc.s(k); qc.h(k); }
				++k;
			}
			for (size_t i = 0; i < numQubits; ++i) {
				for (size_t j = i + 1; j < numQubits; ++j) {
					if (graph.hasEdge(i, j))
						qc.cz(i, j);
				}
			}
			for (size_t i = 0; i < numQubits; ++i) qc.h(i);
			return qc;
		}

		std::string serialize() const {
			std::stringstream result;
			result << "n=" << numQubits << ":";
			for (const auto& gate : singleQubitLayer) {
				if (gate == BinaryCliffordGates::I) result << 0;
				if (gate == BinaryCliffordGates::H) result << 1;
				if (gate == BinaryCliffordGates::S) result << 2;
				if (gate == BinaryCliffordGates::SH) result << 3;
				if (gate == BinaryCliffordGates::HSH) result << 4;
				if (gate == BinaryCliffordGates::HS) result << 5;
			}
			result << ':';
			for (size_t i = 0; i < numQubits; ++i) {
				for (size_t j = i + 1; j < numQubits; ++j) {
					if (graph.hasEdge(i, j)) {
						result << i << "-" << j << ',';
					}
				}
			}
			return result.str();
		}

		class DeserializationException : public std::runtime_error {
		public:
			using std::runtime_error::runtime_error;
		};

		void deserialize(const std::string& input) {
			const auto size = input.size();
			if (input.size() < 3 || !input.starts_with("n=")) throw DeserializationException("Bad input: Should start with \"n=\"");
			std::istringstream stream(input);
			stream.ignore(2);
			int numQubitsIn{};
			stream >> numQubitsIn;
			if (numQubitsIn != numQubits) throw DeserializationException("Number n of qubits does not match this classes template argument");
			if (size - stream.tellg() < 2 + numQubits) throw DeserializationException("Bad input (input to short)");
			stream.ignore(1);
			for (size_t i = 0; i < numQubits; i++) {
				if (stream.eof()) throw DeserializationException("Bad input (input to short)");
				char c{};
				stream >> c;
				switch (c) {
				case '0': singleQubitLayer[i] = BinaryCliffordGates::I; break;
				case '1': singleQubitLayer[i] = BinaryCliffordGates::H; break;
				case '2': singleQubitLayer[i] = BinaryCliffordGates::S; break;
				case '3': singleQubitLayer[i] = BinaryCliffordGates::SH; break;
				case '4': singleQubitLayer[i] = BinaryCliffordGates::HSH; break;
				case '5': singleQubitLayer[i] = BinaryCliffordGates::HS; break;
				default: throw DeserializationException("bad input");
				}
			}
			graph = Graph<numQubits>{};
			stream.ignore(1);
			while (!stream.eof()) {
				if (size - stream.tellg() < 3) return;
				int i1{};
				int i2{};
				char delim{};
				stream >> i1 >> delim >> i2 >> delim;
				if (i1 < 0 || i1 >= numQubits || i2 < 0 || i2 >= numQubits || i1 == i2) throw DeserializationException("bad edge");
				graph.addEdge(i1, i2);
				graph.addEdge(i2, i1);
			}
		}


		//BinaryOperatorSet<numQubits, numQubits> getStabilizer() const {
		//	BinaryOperatorSet<numQubits, numQubits> stabilizer;
		//	for (size_t i = 0; i < numQubits; ++i) {
		//		stabilizer[i] = transformPauli(BinaryPauliOperator<numQubits>::SingleZ(i));
		//	}
		//	return stabilizer;
		//}

	private:
		void transformThroughHadamardLayer(BinaryPauliOperator<numQubits>& op) const {
			const BinaryCliffordGate H{ 0,1,1,0 };
			for (size_t i = 0; i < numQubits; ++i) {
				op.ops[i] = H * op.ops[i];
			}
		}

		void transformThroughCZ(BinaryPauliOperator<numQubits>& op) const {
			for (size_t i = 0; i < numQubits; ++i) {
				for (size_t j = i + 1; j < numQubits; ++j) {
					if (graph.hasEdge(i, j)) {
						MUBTransforms::applyCZ(op, i, j);
					}
				}
			}
		}

		void transformThroughSingleQubitLayer(BinaryPauliOperator<numQubits>& op) const {
			for (size_t i = 0; i < numQubits; ++i) {
				const auto& gate = singleQubitLayer[i];
				auto& o = op.ops[i];
				auto z = 2 * (gate(0, 1) & gate(1, 0)) * (o[0] & o[1]) // phase flip for each gate that contains hadamard if operator is XZ
					+ o[0].toInt() * (gate(0, 0) & gate(1, 0)) * (2 * gate(1, 1).toInt() - 1) // add phase i for X component if gate is HS or S, negative for HS
					+ o[1].toInt() * (gate(0, 1) & gate(1, 1)) * (2 * gate(1, 0).toInt() - 1); // add phase i for Z component if gate is SH or HSH, negative for HSH

				auto v = (2 * gate(0, 1).toInt() - 1);
				auto w = gate(0, 1);
				auto u = o[1].toInt() * (gate(0, 1) & gate(1, 1)) * (2 * gate(0, 1).toInt() - 1);

				op.phase += z;
				op.ops[i] = gate * op.ops[i];
			}
		}
	};


	template<int n, int m>
	bool testHTCircuit(
		const BinaryMatrix<n, n>& adjacencyMatrix,
		const BinaryVector<n>& axx, const BinaryVector<n>& axz,
		const BinaryVector<n>& azx, const BinaryVector<n>& azz,
		const std::array<BinaryPauliOperator<n>, m>& operators
	) {
		BinaryMatrix<n, n> Axx = diag(axx);
		BinaryMatrix<n, n> Axz = diag(axz);
		BinaryMatrix<n, n> Azx = diag(azx);
		BinaryMatrix<n, n> Azz = diag(azz);
		BinaryMatrix<n, m> R;
		BinaryMatrix<n, m> S;
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				R(i, j) = operators[j].ops[i][0];
				S(i, j) = operators[j].ops[i][1];
			}
		}

		const auto lhs = adjacencyMatrix * (Axx * R + Axz * S) - Azx * R - Azz * S;
		return std::ranges::all_of(lhs, [](auto c) {return c == 0; });
	}





	namespace Latex {

		auto toLatex(const BinaryCliffordGate& gate, int spaceBuffer = 4) {
			const auto gateName = toString(gate);
			return "\\gate{" + gateName + "}" + std::string(std::max(0ULL, spaceBuffer - gateName.size()), ' ');
		}

		template<int n>
		auto toLatex(const HTCircuit <n>& circuit) {
			std::vector<std::string> vs(n);
			for (size_t i = 0; i < n; ++i) {
				vs[i] += "\t & " + toLatex(circuit.singleQubitLayer[i]);
			}
			//vs[0] += " \\slice{} ";

			for (size_t i = 0; i < n; ++i) {
				for (size_t j = i + 1; j < n; ++j) {
					if (circuit.adjacencyMatrix(i, j) == 0) continue;

					vs[i] += " & \\ctrl{" + std::to_string(j - i) + "}   ";
					vs[j] += " & \\control{} ";

					using namespace std::ranges::views;
					for (auto k : (iota(0, n) | filter([i, j](auto a) { return a != i && a != j; }))) {
						vs[k] += " & \\qw        ";
					}

				}
			}
			/*	if (!vs[0].ends_with(" \\slice{} "))
					vs[0] += " \\slice{} ";*/

			for (auto& line : vs) {
				line += " & \\gate{H} & \\qw ";
			}
			std::ostringstream out;
			out << "\\begin{quantikz}\n";
			for (const auto& line : vs) {
				out << line << "\\\\ \n";
			}
			out << "\\end{quantikz}\n";

			return out.str();
		}



		template<int n, int m>
		auto toLatex(const HTCircuit <n>& circuit, const std::array<BinaryPauliOperator<n>, m>& ops) {
			std::vector<std::string> vs(n);
			for (size_t i = 0; i < n; ++i) {
				vs[i] += "\t & " + toLatex(circuit.singleQubitLayer[i]);
			}
			//vs[0] += " \\slice{} ";

			for (size_t i = 0; i < n; ++i) {
				for (size_t j = i + 1; j < n; ++j) {
					if (!circuit.graph.hasEdge(i, j)) continue;
					vs[i] += " & \\ctrl{" + std::to_string(j - i) + "}   ";
					vs[j] += " & \\control{} ";

					using namespace std::ranges::views;
					for (auto k : (iota(0, n) | filter([i, j](auto a) { return a != i && a != j; }))) {
						vs[k] += " & \\qw        ";
					}
				}
			}
			/*	if (!vs[0].ends_with(" \\slice{} "))
					vs[0] += " \\slice{} ";*/

			for (auto& line : vs) {
				line += " & \\gate{H} & \\qw ";
			}
			std::ostringstream out;
			out << "\\begin{quantikz}\n";
			out << "\t\\lstick[wires=" << n << ", brackets = none]{";
			for (size_t i = 0; i < m; ++i) {
				out << ops[i].toString(false);
				if (i != m - 1) out << "\\\\";
			}
			out << "} \n";
			out << vs[0] << "\\rstick[wires=" << n << ", brackets = none]{";

			for (size_t i = 0; i < m; ++i) {
				out << circuit.transformPauli(ops[i]).toString(true);
				if (i != m - 1) out << "\\\\";
			}
			out << "} \\\\ \n";
			for (const auto& line : vs) {
				out << line << "\\\\ \n";
			}
			out << "\\end{quantikz}\n";

			return out.str();
		}
	}



	void test3() {
		using namespace BinaryPauli;
		constexpr BinaryMatrix<3, 3> Γ1 = {
			0,1,0,
			1,0,1,
			0,1,0
		};
		constexpr BinaryMatrix<3, 3> Γ2 = {
			0,0,0,
			0,0,0,
			0,0,0
		};
		constexpr BinaryMatrix<3, 3> Γ3 = {
			0,1,0,
			1,0,0,
			0,0,0
		};
		constexpr BinaryMatrix<3, 3> Γ4 = {
			0,0,0,
			0,0,1,
			0,1,0
		};
		constexpr BinaryMatrix<3, 3> Γ5 = {
			0,0,1,
			0,0,0,
			1,0,0
		};
		constexpr auto set1 = parseMubSet<3, 7>("ZII IIZ IZI ZIZ IZZ ZZZ ZZI"); // HHH      Γ2  000111111000
		constexpr auto set2 = parseMubSet<3, 7>("XII IXI IIX XXI IXX XXX XIX"); // III      Γ2  111000000111
		constexpr auto set3 = parseMubSet<3, 7>("YII IXZ IZX YXZ IYY YYY YZX"); // SHH      Γ4  100011111100

		constexpr auto set4 = parseMubSet<3, 7>("XIZ IYI ZIY XYZ ZYY YYX YIX"); // SH S H   Γ5  010101111110
		constexpr auto set5 = parseMubSet<3, 7>("XZI ZXZ IZY YYZ ZYX YXX XIY"); // HSH I SH Γ1  110101001111
		constexpr auto set6 = parseMubSet<3, 7>("YIZ IYZ ZZY YYI ZXX XXY XZX"); // S HS HS  Γ1  111011111100

		constexpr auto set7 = parseMubSet<3, 7>("XZZ ZYZ ZZX YXI IXY XYX YIY"); // S H SH   Γ1  100011111101
		constexpr auto set8 = parseMubSet<3, 7>("YZZ ZYI ZIX XXZ IYX YXY XZY"); // HS HS I  Γ1  111110110001
		constexpr auto set9 = parseMubSet<3, 7>("YZI ZXI IIY XYI ZXY XYY YZY"); // HS H S   Γ3  101110111001


		testHTCircuit(Γ2, { 0,0,0 }, { 1,1,1 }, { 1,1,1 }, { 0,0,0 }, set1); // HHH      Γ4
		testHTCircuit(Γ2, { 1,1,1 }, { 0,0,0 }, { 0,0,0 }, { 1,1,1 }, set2); // III      Γ4
		testHTCircuit(Γ4, { 1,0,0 }, { 0,1,1 }, { 1,1,1 }, { 1,0,0 }, set3); // SHH      Γ4
		testHTCircuit(Γ5, { 0,1,0 }, { 1,0,1 }, { 1,1,1 }, { 1,1,0 }, set4); // 
		testHTCircuit(Γ1, { 1,1,0 }, { 1,0,1 }, { 0,0,1 }, { 1,1,1 }, set5); // HSH I SH Γ1
		testHTCircuit(Γ1, { 1,1,1 }, { 0,1,1 }, { 1,1,1 }, { 1,0,0 }, set6); // S HS HS  Γ1
		testHTCircuit(Γ1, { 1,0,0 }, { 0,1,1 }, { 1,1,1 }, { 1,0,1 }, set7); // S H SH   Γ1
		testHTCircuit(Γ1, { 1,1,1 }, { 1,1,0 }, { 1,1,0 }, { 0,0,1 }, set8); // HS HS I  Γ1
		testHTCircuit(Γ3, { 1,0,1 }, { 1,1,0 }, { 1,1,1 }, { 0,0,1 }, set9); // HS H S   Γ3
	}
}