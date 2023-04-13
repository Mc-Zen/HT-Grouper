#pragma once
#include "gurobi_c++.h"
#include "graph.h"
#include "binary_pauli.h"
#include "symbolic.h"

#include <cassert>
#include <optional>
#include <string_view>
#include <fstream>
#include <chrono>
#include "formatting.h"

namespace Q {
	void qwe() { /**/ }

	/// @brief Find a Local Clifford (if it exists) that rotates a given stabilizer into a given graph state |Γ〉. 
	/// @param graph    Graph that describes the graph state |Γ〉
	/// @param RS       Stabilizer as a list of Pauli operators
	/// @param verbose  If set to true, the generated equations are printed to stdout
	/// @return         If successfull, a list of symplectic 2x2 matrices, corresponding to the 6 single-qubit Clifford gates
	template<int numQubits>
	std::optional<std::array<BinaryCliffordGate, numQubits>> findHTCircuit(const Graph<numQubits>& graph, auto RS, bool verbose = false) {
		using std::cout;
		//constexpr auto numOperatorsPerSet = pow2(n) - 1;
		constexpr auto numOperatorsPerSet = numQubits; // generator suffices
		constexpr auto m = numOperatorsPerSet;
		constexpr auto numEqs = numQubits * numOperatorsPerSet;

		// Declare diagonal symbol matrices for blocks of symplectic matrix. 
		const auto Axx = Math::diag(generateSymbolVector<numQubits>("axx"));
		const auto Axz = Math::diag(generateSymbolVector<numQubits>("axz"));
		const auto Azx = Math::diag(generateSymbolVector<numQubits>("azx"));
		const auto Azz = Math::diag(generateSymbolVector<numQubits>("azz"));

		Math::Matrix<int, numQubits, numOperatorsPerSet> R;
		Math::Matrix<int, numQubits, numOperatorsPerSet> S;
		for (size_t i = 0; i < numQubits; ++i) {
			for (size_t j = 0; j < numOperatorsPerSet; ++j) {
				R(i, j) = RS[j].x(i);
				S(i, j) = RS[j].z(i);
			}
		}

		const auto& adjacencyMatrix = graph.getAdjacencyMatrix();
		const auto lhs = simplified(adjacencyMatrix.cast<int>() * (Axx * R + Axz * S) + Azx * R + Azz * S);

		if (verbose)
			std::cout << R << S << lhs;

		try {
			GRBEnv env{ true };
			env.set(GRB_IntParam_OutputFlag, verbose);
			env.set("LogFile", "mip1.log");
			env.start();

			GRBModel model{ env }; // has variables, constraints and associated attributes
			model.setObjective(GRBLinExpr{ 0 }, GRB_MINIMIZE);

			// Add quadratic constraints
			std::vector<GRBVar> axxVars, axzVars, azxVars, azzVars;
			const std::vector aVars = { &axxVars,&axzVars,&azxVars,&azzVars };
			for (int i = 0; i < numQubits; ++i) {
				axxVars.push_back(model.addVar(0, 1, 0, GRB_BINARY, "axx" + std::to_string(i)));
				axzVars.push_back(model.addVar(0, 1, 0, GRB_BINARY, "axz" + std::to_string(i)));
				azxVars.push_back(model.addVar(0, 1, 0, GRB_BINARY, "azx" + std::to_string(i)));
				azzVars.push_back(model.addVar(0, 1, 0, GRB_BINARY, "azz" + std::to_string(i)));

				model.addQConstr(axxVars[i] * azzVars[i] + axzVars[i] * azxVars[i] == 1, "qc" + std::to_string(i));
			}


			auto toGRBLinExpression = [&](const Term& t) {
				auto getVar = [&](const Variable& t) -> GRBVar& {
					return aVars[((t.name[1] == 'z') * 2) + (t.name[2] == 'z')]->operator[](t.name[3] - '0');
				};
				GRBLinExpr expr;
				const auto simplifiedTerm = simplified(t);
				for (const auto& c : simplifiedTerm.numbers) {
					expr += c.real();
				}
				for (const auto& c : simplifiedTerm.variables) {
					expr += getVar(c);
				}
				assert(simplifiedTerm.products.empty() && "Expression is not linear");
				for (const auto& c : simplifiedTerm.products) {
					assert(c.sums.empty());
					GRBLinExpr p;

					for (const auto& d : c.numbers) {
						p *= d.real();
					}
					for (const auto& d : c.variables) {
						//p *= getVar(d);
					}
					expr += p;
				}
				return expr;
			};

			// dummy variables (used for rhs to ensure that lhs is multiple of 2 (corresponding to 0 in binary field) for each entry
			std::vector<GRBVar> dummyVars;
			for (int i = 0; i < numEqs; ++i) {
				dummyVars.push_back(model.addVar(-1000, 1000, 0.0, GRB_INTEGER, ""));
			}

			int i{};
			for (const auto& row : lhs) {
				model.addConstr(toGRBLinExpression(row) * 0.5 == dummyVars[i], "i" + std::to_string(i));
				++i;
			}


			model.optimize();
			if (model.get(GRB_IntAttr_Status) != 2) { // failed
				return std::nullopt;
			}
			//model.write("model.lp");
			std::array<BinaryCliffordGate, numQubits> singleQubitLayer;
			for (int i = 0; i < numQubits; ++i) {
				singleQubitLayer[i] = { axxVars[i].get(GRB_DoubleAttr_X) ,axzVars[i].get(GRB_DoubleAttr_X) ,azxVars[i].get(GRB_DoubleAttr_X) ,azzVars[i].get(GRB_DoubleAttr_X) };
				if (verbose)
					cout << "U_" << i << "\n" << singleQubitLayer[i];
			}
			return singleQubitLayer;
		}
		catch (const GRBException& e) {
			cout << "Error code = " << e.getErrorCode() << '\n';
			cout << e.getMessage() << '\n';
		}
		catch (...) {
			cout << "Exception during optimization" << '\n';
		}
		return std::nullopt;
	}
	/// @brief Find a Local Clifford (if it exists) that rotates a given stabilizer into a given graph state |Γ〉. 
	/// @param graph    Graph that describes the graph state |Γ〉
	/// @param RS       Stabilizer as a list of Pauli operators
	/// @param verbose  If set to true, the generated equations are printed to stdout
	/// @return         If successfull, a list of symplectic 2x2 matrices, corresponding to the 6 single-qubit Clifford gates
	std::optional<std::vector<BinaryCliffordGate>> findHTCircuit(const Graph<>& graph, auto RS, bool verbose = false) {
		using std::cout;
		//constexpr auto numOperatorsPerSet = pow2(n) - 1;
		int numQubits = graph.numVertices();
		auto numOperatorsPerSet = RS.size(); // generator suffices
		auto numEqs = numQubits * numOperatorsPerSet;

		// Declare diagonal symbol matrices for blocks of symplectic matrix. 
		const auto Axx = Math::diag(generateSymbolVector(numQubits, "axx"));
		const auto Axz = Math::diag(generateSymbolVector(numQubits, "axz"));
		const auto Azx = Math::diag(generateSymbolVector(numQubits, "azx"));
		const auto Azz = Math::diag(generateSymbolVector(numQubits, "azz"));

		auto t0 = std::chrono::high_resolution_clock::now();

		Math::Matrix<int> R(numQubits, numOperatorsPerSet);
		Math::Matrix<int> S(numQubits, numOperatorsPerSet);
		for (size_t i = 0; i < numQubits; ++i) {
			for (size_t j = 0; j < numOperatorsPerSet; ++j) {
				R(i, j) = RS[j].x(i);
				S(i, j) = RS[j].z(i);
			}
		}

		const auto& adjacencyMatrix = graph.getAdjacencyMatrix();
		const auto lhs = simplified(adjacencyMatrix.cast<int>() * (Axx * R + Axz * S) + Azx * R + Azz * S);

		if (verbose)
			std::cout << R << S << lhs;

		try {
			GRBEnv env{ true };
			env.set(GRB_IntParam_OutputFlag, verbose);
			env.set("LogFile", "mip1.log");
			env.start();

			GRBModel model{ env }; // has variables, constraints and associated attributes
			model.setObjective(GRBLinExpr{ 0 }, GRB_MINIMIZE);

			// Add quadratic constraints
			std::vector<GRBVar> axxVars, axzVars, azxVars, azzVars;
			const std::vector aVars = { &axxVars,&axzVars,&azxVars,&azzVars };
			for (int i = 0; i < numQubits; ++i) {
				axxVars.push_back(model.addVar(0, 1, 0, GRB_BINARY, "axx" + std::to_string(i)));
				axzVars.push_back(model.addVar(0, 1, 0, GRB_BINARY, "axz" + std::to_string(i)));
				azxVars.push_back(model.addVar(0, 1, 0, GRB_BINARY, "azx" + std::to_string(i)));
				azzVars.push_back(model.addVar(0, 1, 0, GRB_BINARY, "azz" + std::to_string(i)));

				model.addQConstr(axxVars[i] * azzVars[i] + axzVars[i] * azxVars[i] == 1, "qc" + std::to_string(i));
			}


			auto toGRBLinExpression = [&](const Term& t) {
				auto getVar = [&](const Variable& t) -> GRBVar& {
					return aVars[((t.name[1] == 'z') * 2) + (t.name[2] == 'z')]->operator[](t.name[3] - '0');
				};
				GRBLinExpr expr;
				const auto simplifiedTerm = simplified(t);
				for (const auto& c : simplifiedTerm.numbers) {
					expr += c.real();
				}
				for (const auto& c : simplifiedTerm.variables) {
					expr += getVar(c);
				}
				assert(simplifiedTerm.products.empty() && "Expression is not linear");
				for (const auto& c : simplifiedTerm.products) {
					assert(c.sums.empty());
					GRBLinExpr p;

					for (const auto& d : c.numbers) {
						p *= d.real();
					}
					for (const auto& d : c.variables) {
						//p *= getVar(d);
					}
					expr += p;
				}
				return expr;
			};

			// dummy variables (used for rhs to ensure that lhs is multiple of 2 (corresponding to 0 in binary field) for each entry
			std::vector<GRBVar> dummyVars;
			for (int i = 0; i < numEqs; ++i) {
				dummyVars.push_back(model.addVar(-1000, 1000, 0.0, GRB_INTEGER, ""));
			}

			int i{};
			for (const auto& row : lhs) {
				model.addConstr(toGRBLinExpression(row) * 0.5 == dummyVars[i], "i" + std::to_string(i));
				++i;
			}

			auto t1 = std::chrono::high_resolution_clock::now();

			model.optimize();
			auto t2 = std::chrono::high_resolution_clock::now();
			//println("Time A: {}, Time B {}", t1 - t0, t2 - t1);
			if (model.get(GRB_IntAttr_Status) != 2) { // failed
				return std::nullopt;
			}
			//model.write("model.lp");
			std::vector<BinaryCliffordGate> singleQubitLayer(numQubits);
			for (int i = 0; i < numQubits; ++i) {
				singleQubitLayer[i] = { axxVars[i].get(GRB_DoubleAttr_X) ,axzVars[i].get(GRB_DoubleAttr_X) ,azxVars[i].get(GRB_DoubleAttr_X) ,azzVars[i].get(GRB_DoubleAttr_X) };
				if (verbose)
					cout << "U_" << i << "\n" << singleQubitLayer[i];
			}
			return singleQubitLayer;
		}
		catch (const GRBException& e) {
			cout << "Error code = " << e.getErrorCode() << '\n';
			cout << e.getMessage() << '\n';
		}
		catch (...) {
			cout << "Exception during optimization" << '\n';
		}
		return std::nullopt;
	}


	class HTCircuitFinder {
		GRBEnv env{ true };
		std::unique_ptr<GRBModel> model;

		std::vector<GRBVar> axxVars, axzVars, azxVars, azzVars;
		Math::Matrix<Q::Term> Axx;
		Math::Matrix<Q::Term> Axz;
		Math::Matrix<Q::Term> Azx;
		Math::Matrix<Q::Term> Azz;

		Math::Matrix<int> R;
		Math::Matrix<int> S;

		Math::Matrix<Q::Term> C; // Axx * R + Axz * S
		Math::Matrix<Q::Term> D; // Azx * R + Azz * S

		int numQubits{};
		int numOperatorsPerSet{};

	public:
		HTCircuitFinder(int numQubits, bool verbose=false) : numQubits(numQubits) {
			env.set(GRB_IntParam_OutputFlag, verbose);
			env.set("LogFile", "mip1.log");
			env.start();

			model = std::make_unique<GRBModel>(env);
			model->setObjective(GRBLinExpr{ 0 }, GRB_MINIMIZE);

			// Declare diagonal symbol matrices for blocks of symplectic matrix. 
			Axx = Math::diag(generateSymbolVector(numQubits, "axx"));
			Axz = Math::diag(generateSymbolVector(numQubits, "axz"));
			Azx = Math::diag(generateSymbolVector(numQubits, "azx"));
			Azz = Math::diag(generateSymbolVector(numQubits, "azz"));
			// Add quadratic constraints
			for (int i = 0; i < numQubits; ++i) {
				axxVars.push_back(model->addVar(0, 1, 0, GRB_BINARY, "axx" + std::to_string(i)));
				axzVars.push_back(model->addVar(0, 1, 0, GRB_BINARY, "axz" + std::to_string(i)));
				azxVars.push_back(model->addVar(0, 1, 0, GRB_BINARY, "azx" + std::to_string(i)));
				azzVars.push_back(model->addVar(0, 1, 0, GRB_BINARY, "azz" + std::to_string(i)));

				model->addQConstr(axxVars[i] * azzVars[i] + axzVars[i] * azxVars[i] == 1, "qc" + std::to_string(i));
			}
		}

		HTCircuitFinder(const HTCircuitFinder&) = delete;
		HTCircuitFinder& operator=(const HTCircuitFinder&) = delete;
		HTCircuitFinder(HTCircuitFinder&&) = default;
		HTCircuitFinder& operator=(HTCircuitFinder&&) = default;

		void setOperators(auto&& RS) {
			numOperatorsPerSet = RS.size();
			R.resize(numQubits, numOperatorsPerSet);
			S.resize(numQubits, numOperatorsPerSet);
			for (size_t i = 0; i < numQubits; ++i) {
				for (size_t j = 0; j < numOperatorsPerSet; ++j) {
					R(i, j) = RS[j].x(i);
					S(i, j) = RS[j].z(i);
				}
			}
			C = Axx * R + Axz * S;
			D = Azx * R + Azz * S;
		}

		/// @brief Find a Local Clifford (if it exists) that rotates a given stabilizer into a given graph state |Γ〉. 
		/// @param graph    Graph that describes the graph state |Γ〉
		/// @param RS       Stabilizer as a list of Pauli operators
		/// @param verbose  If set to true, the generated equations are printed to stdout
		/// @return         If successfull, a list of symplectic 2x2 matrices, corresponding to the 6 single-qubit Clifford gates
		std::optional<std::vector<BinaryCliffordGate>> findHTCircuit(const Graph<>& graph, bool verbose = false) {
			using std::cout;
			auto numEqs = numQubits * numOperatorsPerSet;
			std::vector<std::vector<GRBVar>*> aVars = { &axxVars,&axzVars,&azxVars,&azzVars };


			auto t0 = std::chrono::high_resolution_clock::now();

			const auto& adjacencyMatrix = graph.getAdjacencyMatrix();
			const auto lhs = simplified(adjacencyMatrix.cast<int>() * C + D);
			auto t1 = std::chrono::high_resolution_clock::now();

			if (verbose)
				std::cout << R << S << lhs;

			try {


				auto toGRBLinExpression = [&](const Term& t) {
					auto getVar = [&](const Variable& t) -> GRBVar& {
						return aVars[((t.name[1] == 'z') * 2) + (t.name[2] == 'z')]->operator[](t.name[3] - '0');
					};
					GRBLinExpr expr;
					const auto simplifiedTerm = simplified(t);
					for (const auto& c : simplifiedTerm.numbers) {
						expr += c.real();
					}
					for (const auto& c : simplifiedTerm.variables) {
						expr += getVar(c);
					}
					assert(simplifiedTerm.products.empty() && "Expression is not linear");
					for (const auto& c : simplifiedTerm.products) {
						assert(c.sums.empty());
						GRBLinExpr p;

						for (const auto& d : c.numbers) {
							p *= d.real();
						}
						for (const auto& d : c.variables) {
							//p *= getVar(d);
						}
						expr += p;
					}
					return expr;
				};

				//model->feasRelax(GRB_FEASRELAX_QUADRATIC, true, true, true);
				// dummy variables (used for rhs to ensure that lhs is multiple of 2 (corresponding to 0 in binary field) for each entry
				std::vector<GRBVar> dummyVars;
				for (int i = 0; i < numEqs; ++i) {
					dummyVars.push_back(model->addVar(-1000, 1000, 0.0, GRB_INTEGER, ""));
				}

				std::vector<GRBConstr> constraints;
				int i{};
				for (const auto& row : lhs) {
					constraints.push_back(model->addConstr(toGRBLinExpression(row) * 0.5 == dummyVars[i], "i" + std::to_string(i)));
					++i;
				}

				model->optimize();
				auto t2 = std::chrono::high_resolution_clock::now();
				//println("Time A: {}, Time B {}", t1 - t0, t2 - t1);


				for (auto& var : dummyVars) {
					model->remove(var);
				}
				for (auto& constr : constraints) {
					model->remove(constr);
				}
				if (model->get(GRB_IntAttr_Status) != 2) { // failed
					return std::nullopt;
				}
				//model.write("model.lp");
				std::vector<BinaryCliffordGate> singleQubitLayer(numQubits);
				for (int i = 0; i < numQubits; ++i) {
					singleQubitLayer[i] = { axxVars[i].get(GRB_DoubleAttr_X) ,axzVars[i].get(GRB_DoubleAttr_X) ,azxVars[i].get(GRB_DoubleAttr_X) ,azzVars[i].get(GRB_DoubleAttr_X) };
					if (verbose)
						cout << "U_" << i << "\n" << singleQubitLayer[i];
				}
				return singleQubitLayer;
			}
			catch (const GRBException& e) {
				cout << "Error code = " << e.getErrorCode() << '\n';
				cout << e.getMessage() << '\n';
			}
			catch (...) {
				cout << "Exception during optimization" << '\n';
			}
			return std::nullopt;
		}
	};



}

