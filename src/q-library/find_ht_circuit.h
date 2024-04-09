#pragma once
#include "gurobi_c++.h"
#include "graph.h"
#include "binary_pauli.h"

#include <cassert>
#include <optional>
#include <string_view>
#include <fstream>
#include <chrono>
#include "formatting.h"

namespace Q {
	constexpr void qwe() { /**/ }


	class HTCircuitFinder {
		GRBEnv env{ true };
		std::unique_ptr<GRBModel> model;

		std::vector<GRBVar> axxVars, axzVars, azxVars, azzVars;
		// dummy variables (used for rhs to ensure that lhs is multiple of 2 (corresponding to 0 in binary field) for each entry
		std::vector<GRBVar> dummyVars;
		std::vector<GRBQConstr> quadraticConstraints;


		Math::Matrix<int> R;
		Math::Matrix<int> S;


		int numOperatorsPerSet{};

	public:


		HTCircuitFinder(int numQubits, bool verbose = false) {
			env.set(GRB_IntParam_OutputFlag, verbose);
			env.set("LogFile", "mip1.log");
			env.start();

			model = std::make_unique<GRBModel>(env);
			model->setObjective(GRBLinExpr{ 0 }, GRB_MINIMIZE);

			// Declare diagonal symbol matrices for blocks of symplectic matrix. 
		}

		HTCircuitFinder(const HTCircuitFinder&) = delete;
		HTCircuitFinder& operator=(const HTCircuitFinder&) = delete;
		HTCircuitFinder(HTCircuitFinder&&) = default;
		HTCircuitFinder& operator=(HTCircuitFinder&&) = default;





		/// @brief Find a Local Clifford (if it exists) that rotates a given stabilizer into a given graph state |Γ〉. 
		/// @param graph    Graph that describes the graph state |Γ〉
		/// @param RS       Stabilizer as a list of Pauli operators
		/// @param verbose  If set to true, the generated equations are printed to stdout
		/// @return         If successfull, a list of symplectic 2x2 matrices, corresponding to the 6 single-qubit Clifford gates
		std::optional<std::vector<BinaryCliffordGate>> findHTCircuit(
			const Graph<>& graph,
			const std::vector<Pauli>& paulis,
			bool verbose = false
		) {
			auto numQubits = graph.numVertices();
			auto numPaulis = paulis.size();
			auto numEqs = numQubits * numPaulis;
			const auto& gamma = graph.getAdjacencyMatrix();
			updateSize(numQubits, numPaulis);

			std::vector<GRBConstr> constraints;

			for (int i = 0; i < numQubits; ++i) {
				for (int j = 0; j < numPaulis; ++j) {

					GRBLinExpr expr;
					if (paulis[j].x(i)) expr += azxVars[i];
					if (paulis[j].z(i)) expr += azzVars[i];
					for (int k = 0; k < numQubits; ++k) {
						if (gamma(i, k)) {
							if (paulis[j].x(k)) expr += axxVars[k];
							if (paulis[j].z(k)) expr += axzVars[k];
						}
					}
					constraints.push_back(model->addConstr(expr * 0.5 == dummyVars[i * numPaulis + j]));
				}
			}

			return optimize(constraints, verbose);
		}


		/// @brief Find a Local Clifford (if it exists) that rotates a given stabilizer into a given graph state |Γ〉. 
		/// @param graph    Graph that describes the graph state |Γ〉
		/// @param RS       Stabilizer as a list of Pauli operators
		/// @param verbose  If set to true, the generated equations are printed to stdout
		/// @return         If successfull, a list of symplectic 2x2 matrices, corresponding to the 6 single-qubit Clifford gates
		std::optional<std::vector<BinaryCliffordGate>> findHTCircuit(
			const Graph<>& graph,
			const std::vector<Pauli>& paulis,
			const std::vector<int>& qubits,
			bool verbose = false
		) {
			auto numQubits = qubits.size();
			auto numPaulis = paulis.size();
			auto numEqs = numQubits * numPaulis;
			const auto& gamma = graph.getAdjacencyMatrix();
			updateSize(numQubits, numPaulis);

			std::vector<GRBConstr> constraints;

			for (int i = 0; i < numQubits; ++i) {
				for (int j = 0; j < numPaulis; ++j) {

					GRBLinExpr expr;
					if (paulis[j].x(qubits[i])) expr += azxVars[i];
					if (paulis[j].z(qubits[i])) expr += azzVars[i];
					for (int k = 0; k < numQubits; ++k) {
						if (gamma(qubits[i], qubits[k])) {
							if (paulis[j].x(qubits[k])) expr += axxVars[k];
							if (paulis[j].z(qubits[k])) expr += axzVars[k];
						}
					}
					constraints.push_back(model->addConstr(expr * 0.5 == dummyVars[i * numPaulis + j]));
				}
			}

			return optimize(constraints, verbose);
		}


		std::optional<std::vector<BinaryCliffordGate>> optimize(const std::vector<GRBConstr>& constraints, bool verbose) {

			try {
				model->optimize();
				auto t2 = std::chrono::high_resolution_clock::now();
				//println("Time A: {}, Time B {}", t1 - t0, t2 - t1);


				for (auto& constr : constraints) model->remove(constr);

				if (model->get(GRB_IntAttr_Status) != 2) { // failed
					return std::nullopt;
				}
				//model.write("model.lp");
				std::vector<BinaryCliffordGate> singleQubitLayer(numQubits);
				for (int i = 0; i < numQubits; ++i) {
					singleQubitLayer[i] = { axxVars[i].get(GRB_DoubleAttr_X) ,axzVars[i].get(GRB_DoubleAttr_X) ,azxVars[i].get(GRB_DoubleAttr_X) ,azzVars[i].get(GRB_DoubleAttr_X) };
					if (verbose)
						std::cout << "U_" << i << "\n" << singleQubitLayer[i];
				}
				return singleQubitLayer;
			}
			catch (const GRBException& e) {
				std::cout << "Error code = " << e.getErrorCode() << '\n';
				std::cout << e.getMessage() << '\n';
			}
			catch (...) {
				std::cout << "Exception during optimization" << '\n';
			}
			return std::nullopt;
		}


	private:

		int numQubits{};
		void updateSize(int newNumQubits, int numPaulis) {
			auto numEquations = newNumQubits * numPaulis;
			if (dummyVars.size() < numEquations) {
				for (int i = dummyVars.size(); i < numEquations; ++i) {
					dummyVars.push_back(model->addVar(-1000, 1000, 0.0, GRB_INTEGER));
				}
			}

			if (axxVars.size() < newNumQubits) {
				for (int i = axxVars.size(); i < newNumQubits; ++i) {
					axxVars.push_back(model->addVar(0, 1, 0, GRB_BINARY, "axx" + std::to_string(i)));
					axzVars.push_back(model->addVar(0, 1, 0, GRB_BINARY, "axz" + std::to_string(i)));
					azxVars.push_back(model->addVar(0, 1, 0, GRB_BINARY, "azx" + std::to_string(i)));
					azzVars.push_back(model->addVar(0, 1, 0, GRB_BINARY, "azz" + std::to_string(i)));
				}
			}
			if (numQubits != newNumQubits) {
				if (newNumQubits < numQubits) {
					for (int i = newNumQubits; i < numQubits; ++i) {
						model->remove(quadraticConstraints[i]);
					}
					quadraticConstraints.erase(quadraticConstraints.begin() + newNumQubits, quadraticConstraints.end());
				}
				else {
					for (int i = quadraticConstraints.size(); i < newNumQubits; ++i) {
						auto constraint = model->addQConstr(axxVars[i] * azzVars[i] + axzVars[i] * azxVars[i] == 1, "qc" + std::to_string(i));
						quadraticConstraints.push_back(constraint);
					}
				}
				numQubits = newNumQubits;
			}
		}



	};



}

