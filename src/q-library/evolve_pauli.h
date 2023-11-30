
#pragma once

#include <cassert>
#include "quantum_circuit.h"
#include "pauli.h"

namespace Q {

	auto evolvePauli(const Pauli& pauli, const QuantumCircuit& circuit) {
		assert(pauli.numQubits() == circuit.numQubits);

		Pauli result{ pauli };
		for (const auto& gate : circuit.gates) {
			const auto target = gate.target;
			const auto control = gate.control;
			switch (gate.type) {
				using enum QuantumCircuit::GateType;
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

}


