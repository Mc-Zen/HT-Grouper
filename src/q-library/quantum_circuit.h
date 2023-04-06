#pragma once

#include <vector>
#include <array>
#include "binary_pauli.h"
#include "string_utility.h"


namespace Q {


	/// @brief  Class for quantum Clifford circuits. 
	template<int numQubits>
	class QuantumCircuit {
	public:

		enum class GateType {
			I, X, Y, Z,
			H, S, SDG, CX, CZ, SWAP
		};

		struct Gate {
			GateType type;
			int target{};
			int control{};

			int numQubits() const { return (type == GateType::CX || type == GateType::CZ || type == GateType::SWAP) ? 2 : 1; }
			constexpr friend bool operator==(const Gate& g1, const Gate& g2) = default;
		};

		std::vector<Gate> gates;

		void i(int qubit) { gates.emplace_back(GateType::I, qubit); }
		void x(int qubit) { gates.emplace_back(GateType::X, qubit); }
		void y(int qubit) { gates.emplace_back(GateType::Y, qubit); }
		void z(int qubit) { gates.emplace_back(GateType::Z, qubit); }
		void h(int qubit) { gates.emplace_back(GateType::H, qubit); }
		void s(int qubit) { gates.emplace_back(GateType::S, qubit); }
		void sdg(int qubit) { gates.emplace_back(GateType::SDG, qubit); }
		void cx(int control, int target) { gates.emplace_back(GateType::CX, target, control); }
		void cz(int control, int target) { gates.emplace_back(GateType::CZ, target, control); }
		void swap(int qubit1, int qubit2) { gates.emplace_back(GateType::SWAP, qubit1, qubit2); }

		void h(const std::initializer_list<int>& qubits) { for (auto qubit : qubits) gates.emplace_back(GateType::H, qubit); }

		void clear() { gates.clear(); }

		void append(const QuantumCircuit& other) {
			std::ranges::copy(other.gates, std::back_inserter(gates));
		}


		auto transformPauli(const BinaryPauliOperator<numQubits>& input) const {
			BinaryPauliOperator<numQubits> result{ input };
			for (const auto& gate : gates) {
				const auto target = gate.target;
				const auto control = gate.control;
				switch (gate.type) {
					using enum GateType;
				case I: break;
				case X: Clifford::x(result, target); break;
				case Y: Clifford::y(result, target); break;
				case Z: Clifford::z(result, target); break;
				case H: Clifford::h(result, target); break;
				case S: Clifford::s(result, target); break;
				case SDG: Clifford::sdg(result, target); break;
				case CX: Clifford::cx(result, control, target); break;
				case CZ: Clifford::cz(result, control, target); break;
				case SWAP: Clifford::swap(result, control, target); break;
				default:break;
				}

			}
			return result;
		}

		void invert() {
			std::ranges::reverse(gates);
			std::ranges::for_each(gates, [](auto& gate) {
				if (gate.type == GateType::S) gate.type = GateType::SDG;
				else if (gate.type == GateType::SDG) gate.type = GateType::S;
				});
		}
		
		auto inverse() const {
			auto copy = *this;
			copy.invert();
			return copy;
		}

		BinaryOperatorSet<numQubits, numQubits> getStabilizer() const {
			BinaryOperatorSet<numQubits, numQubits> stabilizer;
			for (size_t i = 0; i < numQubits; ++i) {
				stabilizer[i] = transformPauli(BinaryPauliOperator<numQubits>::SingleZ(i));
			}
			return stabilizer;
		}

		std::string serialize() const {
			std::string result;
			for (const auto& gate : gates) {
				const auto target = gate.target;
				const auto control = gate.control;
				switch (gate.type) {
					using enum GateType;
				case I: result += "i(" + std::to_string(target) + ')'; break;
				case X: result += "x(" + std::to_string(target) + ')'; break;
				case Y: result += "y(" + std::to_string(target) + ')'; break;
				case Z: result += "z(" + std::to_string(target) + ')'; break;
				case H: result += "h(" + std::to_string(target) + ')'; break;
				case S: result += "s(" + std::to_string(target) + ')'; break;
				case SDG: result += "sdg(" + std::to_string(target) + ')'; break;
				case CX: result += "cx(" + std::to_string(control) + ',' + std::to_string(target) + ')'; break;
				case CZ: result += "cz(" + std::to_string(control) + ',' + std::to_string(target) + ')'; break;
				case SWAP: result += "swap(" + std::to_string(control) + ',' + std::to_string(target) + ')'; break;
				default:break;
				}
				result += ' ';
			}
			return result;
		}


		struct DeserializationError : public std::runtime_error {
			using runtime_error::runtime_error;
		};


		void deserialize(std::string_view input) {
			clear();
			const auto instructions = split(trim(input), ' ');
			for (const auto& instruction : instructions) {
				if (!instruction.ends_with(')')) throw DeserializationError("Wrong instruction format: missing \")\"");
				const auto openingBracePosition = instruction.find('(');
				if (openingBracePosition == std::string::npos) throw DeserializationError("Wrong instruction format: missing \")\"");
				const auto instructionName = instruction.substr(0, openingBracePosition);
				const auto qubitsStrings = split(instruction.substr(openingBracePosition + 1, instruction.length() - openingBracePosition - 1), ",");
				std::vector<int> qubits;
				std::transform(qubitsStrings.begin(), qubitsStrings.end(), std::back_inserter(qubits), [](const auto& str) {return std::stoi(trim(str)); });
				if (qubits.empty()) throw DeserializationError("Error: no qubit parameter in instruction");
				const int firstQubit = qubits[0];
				if (instructionName == "i") i(firstQubit);
				else if (instructionName == "x") x(firstQubit);
				else if (instructionName == "y") y(firstQubit);
				else if (instructionName == "z") z(firstQubit);
				else if (instructionName == "h") h(firstQubit);
				else if (instructionName == "s") s(firstQubit);
				else if (instructionName == "sdg") sdg(firstQubit);
				else if (instructionName == "cx") {
					if (qubits.size() != 2) throw DeserializationError("The operation cx needs two qubits");
					cx(firstQubit, qubits[1]);
				}
				else if (instructionName == "cz") {
					if (qubits.size() != 2) throw DeserializationError("The operation cx needs two qubits");
					cz(firstQubit, qubits[1]);
				}
				else if (instructionName == "swap") {
					if (qubits.size() != 2) throw DeserializationError("The operation cx needs two qubits");
					swap(firstQubit, qubits[1]);
				}
			}
		}

		friend bool operator==(const QuantumCircuit& qc1, const QuantumCircuit& qc2) = default;



		std::string toString() const {
			auto getGateSymbol = [](const Gate& gate) {
				switch (gate.type) {
					using enum GateType;
				case I: return 'I';
				case X: return 'X';
				case Y: return 'Y';
				case Z: return 'Z';
				case S: return 'S';
				case SDG: return 'D';
				case H: return 'H';
				}
				return '-';
			};

			auto fillUpTo = [](auto& stack, int size) {
				int diff = size - stack.size();
				for (int k = 0; k < diff; ++k) {
					stack.emplace_back('-');
				}
			};

			std::vector<std::vector<char>> matrix(numQubits);
			for (const auto&gate : gates) {
				if (gate.numQubits() == 1) {
					matrix[gate.target].emplace_back(getGateSymbol(gate));
				}
				else {
					auto& stack1 = matrix[gate.target];
					auto& stack2 = matrix[gate.control];
					auto& smallerStack = stack1.size() < stack2.size() ? stack1 : stack2;
					fillUpTo(smallerStack, std::max(stack1.size(), stack2.size()));

					int firstQubit = std::min(gate.target, gate.control);
					int secondQubit = std::max(gate.target, gate.control);
					while (true) {
						for (int i = firstQubit + 1; i < secondQubit; ++i) {
							if (matrix[i].size() >= stack1.size()+1) {
								auto c = matrix[i][stack1.size()];
								if (c != '-') {
									stack1.emplace_back('-');
									stack1.emplace_back('-');
									stack2.emplace_back('-');
									stack2.emplace_back('-');
									continue;
								}
							}
							else {
								fillUpTo(matrix[i], stack1.size() + 1);
							}
						}
						break;
					}
					switch (gate.type) {
						using enum GateType;
					case CX:
						stack1.emplace_back('+');
						stack2.emplace_back('o');
						break;
					case CZ:
						stack1.emplace_back('o');
						stack2.emplace_back('o');
						break;
					case SWAP:
						stack1.emplace_back('x');
						stack2.emplace_back('x');
						break;
					}
				}
			}
			auto maxSize = std::ranges::max_element(matrix, [](const auto& v1, const auto& v2) {return v1.size() < v2.size(); })->size();
			for (auto& stack : matrix) {
				int rest = maxSize - stack.size();
				for (int i = 0; i < rest; ++i) {
					stack.emplace_back('-');
				}
			}
			std::vector<std::string> strings(numQubits);
			std::vector<std::string> strings2(numQubits -1);
			for (int i = 0; i < maxSize; ++i) {
				for (auto& s : strings) s += "--";
				for (auto& s : strings2) s += "  ";
				bool twoGate{};
				for (int j = 0; j < numQubits; ++j) {
					char c = matrix[j][i];
					bool hasTwoGate{};
					if (c == 'o' || c == '+' || c == 'x') {
						hasTwoGate = true;
						twoGate = !twoGate;
					}
					strings[j] += c;
					if (j == numQubits - 1) continue;
					if (twoGate) {
						strings2[j] += '|';
						if(!hasTwoGate)
						strings[j][strings[j].size() - 1] = '|';
					}
					else {
						strings2[j] += ' ';
					}
				}
			}

			for (auto& s : strings) s += "--";
			for (auto& s : strings2) s += "  ";

			std::string result = strings[0] + '\n';
			for (int j = 1; j < numQubits; ++j) {
				result += strings2[j - 1] + '\n';
				result += strings[j] + '\n';
			}
			return result;
		}
	};
}