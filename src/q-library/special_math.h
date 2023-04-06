#pragma once
#include <concepts>
#include <cassert>

namespace Q {
	template<std::integral T>
	constexpr bool isPowerOf2(T x) {
		return !(x & (x - 1)) && x > 0;
	}


	namespace detail {
		struct coeffs
		{
			static constexpr int deBruijnBitPosition2[32] = {
				0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
				31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
			};
		};
	}

	template<std::integral T>
		requires(sizeof(T) == 4) constexpr T log2OfPowerOf2(T x) {
		assert(isPowerOf2(x) && "This function is only designed for numbers which are a power of 2");
		// static const int deBruijnBitPosition2[32] = {
		//	0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
		//	31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
		// };
		return detail::coeffs::deBruijnBitPosition2[(uint32_t)(x * 0x077CB531U) >> 27];
	}



	constexpr uint64_t pow4(uint64_t n) {
		if (n == 0) return 1;
		return (4ULL << 2 * (n - 1));
	}
	

	constexpr uint64_t pow2(uint64_t n) {
		if (n == 0) return 1;
		return 1ULL << n;
	}


	template<std::integral T>
		requires(sizeof(T) == 4) constexpr T bitReverse(T x, int nb) {
		assert(nb > 0 && 32 > nb && "invalid bit count");
		x = (x << 16) | (x >> 16);
		x = ((x & 0x00FF00FF) << 8) | ((x & 0xFF00FF00) >> 8);
		x = ((x & 0x0F0F0F0F) << 4) | ((x & 0xF0F0F0F0) >> 4);
		x = ((x & 0x33333333) << 2) | ((x & 0xCCCCCCCC) >> 2);
		x = ((x & 0x55555555) << 1) | ((x & 0xAAAAAAAA) >> 1);

		return ((x >> (32 - nb)) & (0xFFFFFFFF >> (32 - nb)));
	}



	constexpr int binomialCoeff(int n, int k) {
		int result = 1;

		// C(n, k) = C(n, n-k)
		if (k > n - k) {
			k = n - k;
		}
		for (int i = 0; i < k; ++i) {
			result *= (n - i);
			result /= (i + 1);
		}
		return result;
	}



	template<std::forward_iterator It, std::floating_point T = double>
	std::pair<T, T> meanAndStandardDeviation(It first, It last) {
		using data_type = typename std::iterator_traits<It>::value_type;
		data_type m{}, sd{};
		size_t count{};
		for (auto it = first; it != last; ++it, ++count) {
			m += *it;
		}
		const T mean = m / static_cast<data_type>(count);
		for (auto it = first; it != last; ++it) {
			const auto tmp = *it - mean;
			sd += tmp * tmp;
		}
		return { mean, std::sqrt(sd / static_cast<data_type>(count)) };
	}


}
