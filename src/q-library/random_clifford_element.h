#pragma once

// Implementation from https://aip.scitation.org/doi/abs/10.1063/1.4903507?journalCode=jmp

#include "binary_pauli.h"
#include "special_math.h"
#include "efficient_binary_math.h"
#include <bit>
#include <vector>
#include <ranges>
#include <iostream>

namespace Q {

	uint64_t cliffordGroupSizeModuloPauli(int n) {
		assert(n < 6);
		uint64_t u = 1;
		for (int i = 1; i <= n; ++i) {
			u *= pow4(i) - 1;
		}
		return pow2(n * n) * u;
	}

	template<class T, int n, int m, int p, int q>
	auto directSum(const Math::Matrix<T, m, n>& a, const Math::Matrix<T, p, q>& b) {
		Math::Matrix<T, m + p, n + q> result;
		result.block(0, 0, m, n) = a;
		result.block(m, n, p, q) = b;
		return result;
	}

	template<int n>
	auto createBinaryVector(uint64_t bitstring) {
		BinaryVector<n> output;
		for (size_t i = 0; i < n; ++i) {
			output[i] = bitstring & 1;
			bitstring >>= 1;
		}
		return output;
	}

	template<int n>
	auto toBitstringInteger(const BinaryVector<n>& vec) {
		uint64_t output{};
		for (size_t i = 0; i < n; ++i) {
			output |= vec[i] << i;
		}
		return output;
	}


	template<int n>
	auto transvection(const BinaryVector <n>& k, const BinaryVector<n>& v) {
		return v + k * symplecticInnerProduct(k, v);
	}

	template<int n>
	auto symplecticInnerProduct(const BinaryVector <n>& v, const BinaryVector<n>& w) {
		static_assert(!(n & 1), "n has to be even");
		Binary t = 0;
		for (size_t i = 0; i < n / 2; ++i) {
			t += v[2 * i] * w[2 * i + 1];
			t += w[2 * i] * v[2 * i + 1];
		}
		return t;
	}

	template<int n>
	auto findTransvection(const BinaryVector <n>& x, const BinaryVector<n>& y) {
		static_assert(!(n & 1), "n has to be even");
		BinaryMatrix<n, 2> output;
		if (x == y) return output;
		if (symplecticInnerProduct(x, y) == 1) {
			output.col(0) = x + y;
			return output;
		}
		//    # find a pair where they are both not 00
		BinaryVector<n> z;
		for (size_t i = 0; i < n / 2; ++i) {
			const auto twoI = 2 * i;
			if ((x[twoI].toInt() + x[twoI + 1].toInt()) != 0 && (y[twoI].toInt() + y[twoI + 1].toInt()) != 0) {
				z[twoI] = x[twoI] + y[twoI];
				z[twoI + 1] = x[twoI + 1] + y[twoI + 1];
				if (z[twoI].toInt() + z[twoI + 1].toInt() == 0) {
					z[twoI + 1] = 1;
					if (x[twoI] != x[twoI + 1]) z[twoI] = 1;
				}
				output.col(0) = x + z;
				output.col(1) = y + z;
				return output;
			}
		}
		// didnt find a pair. So look for two places where x has 00 and y doesnt and vice versa
		for (size_t i = 0; i < n / 2; ++i) {
			const auto twoI = 2 * i;
			if ((x[twoI].toInt() + x[twoI + 1].toInt()) != 0 && (y[twoI].toInt() + y[twoI + 1].toInt()) == 0) {
				if (x[twoI] == x[twoI + 1]) {
					z[twoI + 1] = 1;
				}
				else {
					z[twoI + 1] = x[twoI];
					z[twoI] = x[twoI + 1];
				}
				break;
			}
		}
		// didnt find a pair. So look for two places where x has 00 and y doesnt and vice versa
		for (size_t i = 0; i < n / 2; ++i) {
			const auto twoI = 2 * i;
			if ((x[twoI].toInt() + x[twoI + 1].toInt()) == 0 && (y[twoI].toInt() + y[twoI + 1].toInt()) != 0) {
				if (y[twoI] == y[twoI + 1]) {
					z[twoI + 1] = 1;
				}
				else {
					z[twoI + 1] = y[twoI];
					z[twoI] = y[twoI + 1];
				}
				break;
			}
		}
		output.col(0) = x + z;
		output.col(1) = y + z;
		return output;
	}

	template<int n>
	auto symplectic(uint64_t i) {
		constexpr auto twoN = 2 * n;
		const auto s = ((1 << twoN) - 1);
		const auto k = (i % s) + 1;
		i /= s;

		using Vec = BinaryVector<twoN>;
		Vec f1 = createBinaryVector<twoN>(k);
		Vec e1{ 1 };
		const auto T = findTransvection(e1, f1);
		const auto bits = createBinaryVector<twoN - 1>(i % (1 << (twoN - 1)));

		auto eprime = e1;
		for (size_t i = 2; i < twoN; ++i) eprime[i] = bits[i - 1];
		//std::cout << T<<eprime;
		auto h0 = transvection<twoN>(T.col(0), eprime);
		h0 = transvection<twoN>(T.col(1), h0);

		if (bits[0] == 1) {
			f1 *= 0;
		}

		const auto id2 = BinaryMatrix<2, 2>::identity();
		BinaryMatrix<twoN, twoN> g;
		if constexpr (n == 1) {
			g = id2;
		}
		else {
			g = directSum(id2, symplectic<n - 1>(i >> (twoN - 1)));
		}

		for (size_t j = 0; j < twoN; ++j) {
			g.col(j) = transvection(Vec{ T.col(0) }, Vec{ g.col(j) });
			g.col(j) = transvection(Vec{ T.col(1) }, Vec{ g.col(j) });
			g.col(j) = transvection(Vec{ h0 }, Vec{ g.col(j) });
			g.col(j) = transvection(Vec{ f1 }, Vec{ g.col(j) });
		}
		return g;
	}








	/// @brief Much faster implementation (about 20x) by using a lot of bit hacks
	namespace efficient {
		template<int n>
		using BVector = uint_fast16_t;

		template<int m, int n>
		struct BMatrix {
			static_assert(m < 16);
			std::array<BVector<m>, n> cols{};
		};

		template<int m, int n>
		auto toMatrix(const BMatrix<m, n>& mat) {
			BinaryMatrix<m, n> result;
			for (size_t i = 0; i < n; ++i) {
				result.col(i) = createBinaryVector<m>(mat.cols[i]);
			}
			return result;
		}

		template<int n>
		int symplecticInnerProduct(BinaryVector<n> v, BinaryVector<n> w) {
			static_assert(!(n & 1), "n has to be even");
			constexpr uint64_t mask = 0b1010101010101010101010101010101010101010101010101010101010101010;
			return (std::popcount(((v() << 1) & w()) & mask) + std::popcount(((w() << 1) & v()) & mask)) & 1;
		}

		template<int n>
		auto transvection(BinaryVector<n> k, BinaryVector<n> v) {
			return v + (k() * symplecticInnerProduct<n>(k, v));
		}


		template<int n>
		auto findTransvection(BinaryVector<n> x, BinaryVector<n> y) {
			BinaryColMatrix<n, 2> output;
			if (x == y) return output;
			if (symplecticInnerProduct<n>(x, y) == 1) {
				output.cols[0] = x + y;
				return output;
			}


			//    # find a pair where they are both not 00
			constexpr uint64_t mask = 0b1010101010101010101010101010101010101010101010101010101010101010;
			const auto l = (x() | (x() << 1)) & (y() | (y() << 1)) & mask;
			const auto iii = std::countr_zero(l) - 1;

			if (iii != 8 * sizeof(l) - 1) {
				uint64_t z{};
				const auto xx = x() >> iii;
				const auto yy = y() >> iii;
				z = (xx ^ yy) & 3;
				if (!z) {
					z = 2 | static_cast<uint64_t>((xx & 1) != (xx & 2)); // ??
				}
				z <<= iii;
				output.cols[0] = x + z;
				output.cols[1] = y + z;
				return output;
			}

			//for (size_t i = 0; i < n / 2; ++i) {
			//	const auto twoI = 2 * i;
			//	const auto xx = x >> twoI;
			//	const auto yy = y >> twoI;
			//	if (xx & 3 && yy & 3) {
			//		BVector<n> z{};
			//		z = (xx ^ yy) & 3;
			//		if (!z) {
			//			z = 2 | static_cast<uint_fast16_t>((xx & 1) != (xx & 2)); // ??
			//		}
			//		//if (twoI != iii) std::terminate();
			//		//assert(twoI == iii);
			//		z <<= twoI;
			//		output.cols[0] = x ^ z;
			//		output.cols[1] = y ^ z;
			//		return output;
			//	}
			//}
			//if (iii != 8 * sizeof(iii) - 1) std::terminate();
			uint64_t z1{};
			// didnt find a pair. So look for two places where x has 00 and y doesnt and vice versa
			for (size_t i = 0; i < n / 2; ++i) {
				const auto twoI = 2 * i;
				const auto xx = x() >> twoI;
				const auto yy = y() >> twoI;
				if ((xx & 3) && !(yy & 3)) {
					if ((xx & 1) == ((xx & 2) >> 1)) {
						z1 = 2;
					}
					else {
						z1 = ((xx & 1) << 1) | ((xx & 2) >> 1);
					}
					z1 <<= twoI;
					break;
				}
			}
			uint64_t z2{};
			// didnt find a pair. So look for two places where x has 00 and y doesnt and vice versa
			for (size_t i = 0; i < n / 2; ++i) {
				const auto twoI = 2 * i;
				const auto xx = x() >> twoI;
				const auto yy = y() >> twoI;
				if (!(xx & 3) && (yy & 3)) {
					if ((yy & 1) == ((yy & 2) >> 1)) {
						z2 = 2;
					}
					else {
						z2 = ((yy & 1) << 1) | ((yy & 2) >> 1);
					}
					z2 <<= twoI;
					break;
				}
			}
			const auto z = z1 | z2;
			output.cols[0] = x + z;
			output.cols[1] = y + z;
			return output;
		}


		template<int n>
		auto symplectic(uint64_t i) {
			constexpr auto twoN = 2 * n;
			const auto s = ((1 << twoN) - 1);
			const auto k = (i % s) + 1;
			i /= s;

			using Vec = BinaryVector<twoN>;
			Vec f1 = k;
			const Vec e1 = 1;
			const auto T = findTransvection<twoN>(e1, f1);
			const auto bits = i % (1 << (twoN - 1));

			const auto eprime = ((bits << 1) & ~3) | 1;
			//print("{:0>8b}, {:0>8b}, {:0>8b}\n", T.cols[0], T.cols[1], eprime);

			const auto h0 = transvection<twoN>(T.cols[1], transvection<twoN>(T.cols[0], eprime));

			if (bits & 1) {
				f1 *= 0;
			}

			BinaryColMatrix<twoN, twoN> g;
			if constexpr (n == 1) {
				g.cols[0] = 1; // set to 2x2 identity matrix
				g.cols[1] = 2;
			}
			else {
				const auto m = symplectic<n - 1>(i >> (twoN - 1));
				g.cols[0] = 1;
				g.cols[1] = 2;
				for (int j = 0; j < 2 * n - 2; ++j) {
					g.cols[j + 2] = m.cols[j]() << 2;
				}
			}

			for (size_t j = 0; j < twoN; ++j) {
				g.cols[j] = transvection<twoN>(T.cols[0], g.cols[j]);
				g.cols[j] = transvection<twoN>(T.cols[1], g.cols[j]);
				g.cols[j] = transvection<twoN>(h0, g.cols[j]);
				g.cols[j] = transvection<twoN>(f1, g.cols[j]);
			}
			return g;
		}



		template<int numQubits>
		using RowClifford = BinaryRowMatrix<numQubits, numQubits>;
		template<int numQubits>
		using ColClifford = BinaryColMatrix<numQubits, numQubits>;

		template<int numQubits>
		struct Clifford {
			BinaryRowMatrix<numQubits, numQubits> Axx;
			BinaryRowMatrix<numQubits, numQubits> Axz;
			BinaryRowMatrix<numQubits, numQubits> Azx;
			BinaryRowMatrix<numQubits, numQubits> Azz;
		};


		/// @brief Generate a clifford object of 4 (n x n) block matrices from a (2n x 2n) symplectic matrix. 
		///        NOTE: This version does a transposition (the result is still a clifford operation). 
		template<int numQubitsTimes2>
		auto cliffordFrom2n2nSymplectic(const BinaryColMatrix<numQubitsTimes2, numQubitsTimes2>& symplecticMatrix) {
			constexpr auto numQubits = numQubitsTimes2 / 2;
			Clifford<numQubits> clifford;
			constexpr uint64_t maskX = (1 << numQubits) - 1; // first n bits
			constexpr uint64_t maskZ = maskX << numQubits; // second n bits


			for (int row = 0; row < numQubits; ++row) {
				for (int col = 0; col < numQubits; ++col) {
					clifford.Axx[col].set(row, symplecticMatrix[2 * col].get(2 * row));
					clifford.Axz[col].set(row, symplecticMatrix[2 * col].get(2 * row + 1));
					clifford.Azx[col].set(row, symplecticMatrix[2 * col + 1].get(2 * row));
					clifford.Azz[col].set(row, symplecticMatrix[2 * col + 1].get(2 * row + 1));
				}
			}
			//for (int i = 0; i < numQubits; ++i) {
			//	clifford.Axx[i] = maskX & symplecticMatrix[i]();
			//	clifford.Axz[i] = (maskZ & symplecticMatrix[i]()) >> numQubits;
			//	clifford.Azx[i] = (maskX & symplecticMatrix[i + numQubits]());
			//	clifford.Azz[i] = (maskZ & symplecticMatrix[i + numQubits]()) >> numQubits;
			//}
			//println("{}", symplecticMatrix);
			//for (int i = 0; i < 2*numQubits; ++i) {
			//	println("{}", symplecticMatrix[i]);
			//}
			//println("");
			//println("");
			//for (int i = 0; i < numQubits; ++i) {
			//	println("{} {}", clifford.Axx[i], clifford.Axz[i]);
			//}
			//println("");
			//for (int i = 0; i < numQubits; ++i) {
			//	println("{} {}", clifford.Azx[i], clifford.Azz[i]);
			//}
			//println("");
			//println("");
			//println("");

			return clifford;
		}

	}

}