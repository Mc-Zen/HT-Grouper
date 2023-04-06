
#pragma once


#include "get_H_base_from_irreducible_polynomial.h"
#include "binary_pauli.h"
#include "mub.h"


namespace Q {

	template<int n>
	auto generateMub(const std::vector<uint64_t>& base) {
		constexpr auto d = pow2(n);
		assert(base.size() == d && "The base needs d-1 elements with d=2^n");
		BinaryOperatorSet<n, d - 1> xOperators;
		BinaryOperatorSet<n, d - 1> zOperators;


		for (size_t q = 0; q < d - 1; ++q) {
			BinaryPauliOperator<n> op;
			const auto element = base[q + 1];
			for (size_t j = 0; j < n; ++j) {
				if (element & (1ULL << j)) {
					op[j] = BinaryPauli::X;
				}
			}
			xOperators[q] = op;
			//print("{} {}\n", op.toString(), element);
		}

		auto tr = [](auto i) {
			return i & (1ULL << (n - 1));
		};

		for (size_t q = 0; q < d - 1; ++q) {
			BinaryPauliOperator<n> op;
			for (size_t k = 0; k < n; ++k) {
				auto qpk = (q + k+d-2) % (d - 1) + 1;
				if (tr(base[qpk])) {
					op[k] = BinaryPauli::Z;
				}
			}
			zOperators[q] = op;
			//print("{} {}\n", op.toString(), element);
		}

		Mub<n> mub;
		mub.push_back(zOperators);
		mub.push_back(xOperators);
		for (size_t i = 0; i < d - 1; ++i) {
			MubSet<n> mubSet;
			for (size_t j = 0; j < d - 1; ++j) {
				BinaryPauliOperator<n> op = xOperators[j];
				op *= zOperators[(j + i) % (d-1)];
				op.resetPhaseToTreatXZasY();
				mubSet[j] = op;
			}
			mub.push_back(mubSet);
		}
		//mub[0][2] = { "XXYX" };
		//printMUB(mub);

		//std::cout << "Is commuting: " << areMubwiseCommuting(mub) <<", is each operator unique: ";
		//std::cout << isMub(mub) << "\n";
		return std::move(mub);
	}

}
