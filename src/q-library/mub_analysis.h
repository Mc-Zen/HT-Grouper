#pragma once

#include "binary_pauli.h"
#include "mub.h"
#include <vector>
#include <variant>
#include <format>

namespace Q {

	/// Specific for 2-qubit MUBs
	namespace Mub2 {

		auto getMub() {
			std::vector mub = {
				parseMubSet<2, 3>("XY YZ ZX"),
				parseMubSet<2, 3>("XX YY ZZ"),
				parseMubSet<2, 3>("XI XZ IZ"),
				parseMubSet<2, 3>("YI YX IX"),
				parseMubSet<2, 3>("ZI ZY IY")
			};
			return mub;
		}
	}


	/// Specific for 3-qubit MUBs
	namespace Mub3 {

		auto get234Mub() {
			const std::vector mub = {
				parseMubSet<3, 7>("ZII IIZ IZI ZIZ IZZ ZZZ ZZI") ,
				parseMubSet<3, 7>("XII IXI IIX XXI IXX XXX XIX"),
				parseMubSet<3, 7>("YII IXZ IZX YXZ IYY YYY YZX"),
				parseMubSet<3, 7>("XIZ IYI ZIY XYZ ZYY YYX YIX"),
				parseMubSet<3, 7>("XZI ZXZ IZY YYZ ZYX YXX XIY"),
				parseMubSet<3, 7>("YIZ IYZ ZZY YYI ZXX XXY XZX"),
				parseMubSet<3, 7>("XZZ ZYZ ZZX YXI IXY XYX YIY"),
				parseMubSet<3, 7>("YZZ ZYI ZIX XXZ IYX YXY XZY"),
				parseMubSet<3, 7>("YZI ZXI IIY XYI ZXY XYY YZY")
			};
			return mub;
		}

		auto get090Mub() {
			auto mub = get234Mub();
			std::ranges::for_each(mub, [](auto& base) {
				MUBTransforms::localPermutationXYZ(base, 0);
				MUBTransforms::localYZSwap(base, 1);
				MUBTransforms::localYZSwap(base, 2);
				MUBTransforms::applyCZ(base, 1, 2);
				}
			);
			return mub;
		}

		auto get162Mub() {
			auto mub = get234Mub();
			std::ranges::for_each(mub, [](auto& base) {
				MUBTransforms::localPermutationXYZ(base, 0);
				MUBTransforms::localYZSwap(base, 1);
				MUBTransforms::localYZSwap(base, 2);
				MUBTransforms::applyCZ(base, 0, 1);
				}
			);
			return mub;
		}

		auto get306Mub() {
			auto mub = get162Mub();
			std::ranges::for_each(mub, [](auto& base) {
				MUBTransforms::localYZSwap(base, 0);
				MUBTransforms::localYZSwap(base, 1);
				MUBTransforms::localXYSwap(base, 2);
				MUBTransforms::applyCZ(base, 1, 2);
				}
			);
			return mub;
		}
	}


	/// Specific for 4-qubit MUBs
	namespace Mub4 {


		auto get20429Mub() {
			std::vector mub20429 = {
				parseMubSet<4, 15>("ZZZZ ZZIZ ZZZI ZIII                                                       "),
				parseMubSet<4, 15>("XXXX XIII IXII XXIX                                                       "),
				parseMubSet<4, 15>("YYYY YZIZ ZYZI YXIX                                                       "),
				parseMubSet<4, 15>("YYXY YZZI ZXII XYZY                                                       "),
				parseMubSet<4, 15>("YYYX YIII IYZZ YXZX                                                       "),
				parseMubSet<4, 15>("YXXX XZZZ ZXZI XYIX                                                       "),
				parseMubSet<4, 15>("XYYY YIZI IYII YYIX                                                       "),
				parseMubSet<4, 15>("YXYX XZII ZYII YXZY                                                       "),
				parseMubSet<4, 15>("XYXX YZII ZXZZ XXIY                                                       "),
				parseMubSet<4, 15>("YYXX YIZZ IXIZ XYIY                                                       "),
				parseMubSet<4, 15>("YXYY XIIZ IYIZ YXIY                                                       "),
				parseMubSet<4, 15>("XXXY XZIZ ZXIZ XXZX                                                       "),
				parseMubSet<4, 15>("XYXY YIIZ IXZI XXZY                                                       "),
				parseMubSet<4, 15>("YXXY XIZI IXZZ XYZX                                                       "),
				parseMubSet<4, 15>("XXYX XIZZ IYZI YYZY                                                       "),
				parseMubSet<4, 15>("XXYY XZZI ZYZZ YY1Y                                                       "),
				parseMubSet<4, 15>("XYYX YZZZ ZY1Z YYZX                                                       ")
			};
			// Complete MUB as O[r,c] = O[r,c-4] * O[r,c-1]
			for (auto& base : mub20429) {
				for (int64_t i = 4; i < 15; ++i) {
					base[i] = base[i - 4];
					base[i] *= base[i - 1];
					base[i].resetPhaseToTreatXZasY();
				}
			}
			return mub20429;
		}


		// LC-equivalence classes
		enum class GraphType {
			fsep,     // fully separable (graph with no edges), 1+1+1+1
			pair,     // one pair, two isolated vertices, 1+1+2
			triangle, // one isolated vertex, the others are connected (arbitrarily), 1+3
			twoPairs, // two pairs with two vertices each, 2+2
			line,     // connected line graph
			star,     // star graph
		};

		struct PairSpec {
			int isolatedVertex1{}; // holds the lower of the two
			int isolatedVertex2{};
			int numSwaps() const {
				if (isolatedVertex1 == 1 && isolatedVertex2 == 2) return 2;
				if (isolatedVertex2 - isolatedVertex1 == 2) return 1;
				return 0;
			}
		};

		struct TriangleSpec {
			int isolatedVertex{};
		};

		struct TwoPairsSpec {
			int firstPairVertex1; // always holds 0
			int firstPairVertex2; // holds 1, 2 or 3
		};

		struct LineSpec {
			enum class Type {
				u, // 1,2 form an edge and 3,4
				c, // 2,3 form an edge and 1,4
				x  // 1,3 form an edge and 2,4
			};
			Type type{};
		};



		/// @brief Sector length distribution where the first number represents the number of 
		///        operators with Pauli weight 4, the second the number of operators with Pauli
		///        weight 3 and so on. 
		using SLD = std::array<int, 5>;

		/// @brief MUB SLD structure (careful: this is not really the entanglement structure
		///        because the star graph and the graph with two pairs have the same SLD. 
		///        The figure represent the numbers of occurences of the unique SLDs in a MUB:
		///         
		///         index    SLD           Graph                         Separability       #LC-classes
		///         -----------------------------------------------------------------------------------
		///          [0]   1 4 6 4 1    (graph with no edges, 1+1+1+1)        4                 1
		///          [1]   1 2 4 6 3    (one pair, 1+1+2)                     3                 6 
		///          [2]   1 1 3 7 4    (one triple, 1+3)                     2                 4  
		///          [3]   1 0 6 0 9    (two pairs 2+2 or star)               2 / 1             2 / 1
		///          [4]   1 0 2 8 5    (line graph)                          1                 3 (u,c and x)
		using SLDStructure = std::array<int, 5>;



		/// @brief Details about one measurement basis (one stabilizer) of a MUB
		struct SetCharacterisation {
			GraphType graphType;
			SLD sld{};
			std::variant<PairSpec, TwoPairsSpec, LineSpec, TriangleSpec> spec;

			int numCZ() const {
				switch (graphType) {
					using enum Q::Mub4::GraphType;
				case fsep: return 0;
				case pair: return 1;
				case triangle: {
					int isolatedVertex = std::get<TriangleSpec>(spec).isolatedVertex;
					return (isolatedVertex == 1 || isolatedVertex == 2) ? 4 : 2; }
				case twoPairs: return 2;
				case line: return std::get<LineSpec>(spec).type == LineSpec::Type::u ? 4 : std::get<LineSpec>(spec).type == LineSpec::Type::x ? 5 : 3;
				case star: return 3;
				default: throw std::exception{};
				}
			}

			int numSwaps() const {
				switch (graphType) {
					using enum Q::Mub4::GraphType;
				case fsep: return 0;
				case pair: {
					auto s = std::get<PairSpec>(spec);
					return s.numSwaps();
				}
				case triangle: {
					int isolatedVertex = std::get<TriangleSpec>(spec).isolatedVertex;
					return 0;
					//return (isolatedVertex == 1 || isolatedVertex == 2) ? 1 : 0;
				}
				case twoPairs: {
					int firstPairVertex2 = std::get<TwoPairsSpec>(spec).firstPairVertex2;
					if (firstPairVertex2 == 1) return 0;
					if (firstPairVertex2 == 2) return 1;
					return 2;
				}
				case line: {
					auto type = std::get<LineSpec>(spec).type;
					return 0;
					return type == LineSpec::Type::x ? 1 : 0;
				}
				case star: return 0;
				default: throw std::exception{};
				}
			}
		};

		struct MubCharacterisation {
			std::array<SetCharacterisation, 17> setCharacterisations;
			SLDStructure sldStructure{};


			int fullSeparableCount{};
			std::vector<PairSpec> pairs;
			std::vector<TriangleSpec> triangles;
			std::vector<TwoPairsSpec> twoPairs;
			int starCount{};
			int cCount{};
			int uCount{};
			int xCount{};

			int totalNumCZ{};
			int totalNumSwaps{};
			int maxNumCZPerCircuit{};
			int maxNumSwapsPerCircuit{};
			int maxNum2QubitPerCircuit{};

			std::string toString() const {
				auto str = std::format("{}  {}:: {}* {}c {}u {}x ", sldStructure, fullSeparableCount, starCount, cCount, uCount, xCount);
				for (const auto& pair : pairs) {
					str += ":|" + std::to_string(pair.numSwaps()) + ' ';
				}
				for (const auto& triangle : triangles) {
					auto i = triangle.isolatedVertex;
					switch (i) {
					case 0: str += "._ "; break;
					case 1: str += "n_ "; break;
					case 2: str += "_n "; break;
					case 3: str += "_. "; break;
					default: str += "? ";
					}
				}
				for (const auto& tp : twoPairs) {
					switch (tp.firstPairVertex2) {
					case 1: str += "-- "; break;
					case 2: str += "-| "; break;
					case 3: str += "|| "; break;
					default: str += "? ";
					}
				}
				return str + std::format(" #CZ={} #Swap={}, #max-CZ={} #max-Swap={}", totalNumCZ, totalNumSwaps, maxNumCZPerCircuit, maxNumSwapsPerCircuit);
			}

			auto num2QubitGates() const {
				return totalNumSwaps * 3 + totalNumCZ;
			}
		};





		auto graphCharacterizeMub(const Mub<4>& mub) {
			MubCharacterisation characterisation;
			constexpr size_t n = 4;
			constexpr size_t numSets = pow2(n) + 1;


			auto findIsolatedQubit = [](const auto& set) {
				int isolatedQubit{ -1 };
				for (int i = 0; i < 4; ++i) {
					if (countIdentitiesInMubSet(set, i) == 7) isolatedQubit = i;
				}
				if (isolatedQubit == -1) throw std::runtime_error("Illformed MUB");
				return isolatedQubit;
			};
			auto findTwoIsolatedQubits = [](const auto& set) {
				int isolatedQubit1{ -1 };
				int isolatedQubit2{ -1 };
				for (int i = 0; i < 4; ++i) {
					if (countIdentitiesInMubSet(set, i) == 7) {
						if (isolatedQubit1 == -1) isolatedQubit1 = i;
						else isolatedQubit2 = i;
					}
				}
				if (isolatedQubit1 == -1 || isolatedQubit2 == -1) throw std::runtime_error("Illformed MUB");
				return std::make_pair(isolatedQubit1, isolatedQubit2);
			};

			for (size_t i = 0; i < numSets; ++i) {
				auto& setCharacterisation = characterisation.setCharacterisations[i];
				const auto& set = mub[i];
				setCharacterisation.sld = getSLD(set);

				// The first number of the sld (pauli weight 4) uniquely identifies the SLD (for 4 qubits)
				switch (setCharacterisation.sld[0]) {
				case 1:
					++characterisation.sldStructure[0];
					setCharacterisation.graphType = GraphType::fsep;
					++characterisation.fullSeparableCount;
					break;
				case 3: {
					++characterisation.sldStructure[1];
					setCharacterisation.graphType = GraphType::pair;
					const auto [v1, v2] = findTwoIsolatedQubits(set);
					setCharacterisation.spec = PairSpec{ v1, v2 };
					characterisation.pairs.emplace_back(v1, v2);
					break;
				}
				case 4:
					++characterisation.sldStructure[2];
					setCharacterisation.graphType = GraphType::triangle;
					setCharacterisation.spec = TriangleSpec{ findIsolatedQubit(set) };
					characterisation.triangles.emplace_back(findIsolatedQubit(set));
					break;
				case 9: {
					++characterisation.sldStructure[3];
					bool IIAA{}; // does an operator starting with II exist
					bool AIIA{}; // does an operator with II in the middle exist
					for (const auto& op : set) {
						if (op[0] == BinaryPauli::I && op[1] == BinaryPauli::I) IIAA = true;
						if (op[1] == BinaryPauli::I && op[2] == BinaryPauli::I) AIIA = true;
					}
					if (IIAA && AIIA) { // star graphs always all combinations (IIAA, IAIA, IAAI, AIIA, AIAI, AAII)

						setCharacterisation.graphType = GraphType::star;
						++characterisation.starCount;
					}
					else { // two-pair graphs only have (IIAA, AAII) or (IAAI, AIIA) or (IAIA, AIAI)
						setCharacterisation.graphType = GraphType::twoPairs;
						if (IIAA) setCharacterisation.spec = TwoPairsSpec{ 0,1 };
						else if (AIIA) setCharacterisation.spec = TwoPairsSpec{ 0,3 };
						else setCharacterisation.spec = TwoPairsSpec{ 0,2 };
						characterisation.twoPairs.emplace_back(std::get<TwoPairsSpec>(setCharacterisation.spec));
					}
					break;
				}
				case 5: {
					++characterisation.sldStructure[4];
					setCharacterisation.graphType = GraphType::line;
					bool IIAA{}, AIIA{}, IAIA{};
					for (const auto& op : set) {
						if (op[0] == BinaryPauli::I && op[1] == BinaryPauli::I) IIAA = true;
						if (op[1] == BinaryPauli::I && op[2] == BinaryPauli::I) AIIA = true;
						if (op[0] == BinaryPauli::I && op[2] == BinaryPauli::I) IAIA = true;
					}
					if (IIAA) {
						setCharacterisation.spec = LineSpec{ LineSpec::Type::c };
						++characterisation.cCount;
					}
					else if (AIIA) {
						setCharacterisation.spec = LineSpec{ LineSpec::Type::u };
						++characterisation.uCount;
					}
					else if (IAIA) {
						setCharacterisation.spec = LineSpec{ LineSpec::Type::x };
						++characterisation.xCount;
					}
					break;
				}
				default:
					throw std::runtime_error("Invalid SLD encountered");
				}

				const auto numCZ = setCharacterisation.numCZ();
				const auto numSwaps = setCharacterisation.numSwaps();
				characterisation.totalNumCZ += numCZ;
				characterisation.totalNumSwaps += numSwaps;
				characterisation.maxNumCZPerCircuit = std::max(characterisation.maxNumCZPerCircuit, numCZ);
				characterisation.maxNumSwapsPerCircuit = std::max(characterisation.maxNumSwapsPerCircuit, numSwaps);
				characterisation.maxNum2QubitPerCircuit = std::max(characterisation.maxNum2QubitPerCircuit, numSwaps*3+numCZ);
			}

			return characterisation;
		}
	}
}
