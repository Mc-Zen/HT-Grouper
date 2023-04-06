#pragma once

#include <vector>
#include <cassert>
#include <cstdint>
#include "formatting.h"
#include "special_math.h"

namespace Q {

	/// @brief Get a basis of the Hilbert space with dimension 2^n from an irreducible polynomial of degree n. If the 
	///        given polynomial is reducible, an exception is thrown. 
	/// @param degree Degree of the irreducible polynomial. 
	/// @param coeffs Coefficient for the irreducible polynomial through bn*a^n + ... + b2*a^2 + b1*a + b0 so that bn 
	///               is the first entry of coeffs and b0 the last. 
	/// @return A vector with 2^degree bitstrings representing a basis of the Hilbert space. The least significant bit 
	///         represents the coefficient for a^0=1 while the nth bit represents the coefficient for a^n. 
	auto get_H_base_from_irreducible_polynomial(int degree, const std::vector<int>& coeffs) {
		assert(degree > 0 && degree < 64 && "Degree has to be at least one and at max 64");
		assert(coeffs.size() == degree + 1);
		const uint64_t mask = (1ULL << degree) - 1; // Mask with the first (degree-1) bits set to 1

		// when an "overflow" happens, i.e.the nth bit is 1 we want to replace the a^n term according to the given polynomial
		uint64_t overflowMask{}; 
		for (size_t i = 0; i < degree + 1; ++i) {
			overflowMask <<= 1;
			overflowMask |= (coeffs[i] & 1);
		}

		std::vector<uint64_t> base{ 0,1 };
		uint64_t currentElement = 2;
		while (currentElement != 1) {
			base.push_back(currentElement);
			if (base.size() > 1 + pow2(degree)) throw std::runtime_error("Error: polynomial is not irreducible");
			currentElement <<= 1;
			if (currentElement & (1ULL << degree)) {
				currentElement ^= overflowMask;
				currentElement &= mask;
			}
		} 

		if (base.size() != pow2(degree)) throw std::runtime_error("Error: polynomial is not irreducible");
		return base;
	}
}
