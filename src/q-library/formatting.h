
#pragma once
#include "binary.h"
#include "binary_phase.h"
#include <type_traits>
#include <complex>
#include <iostream>
#include <numeric>
#include <vector>
#include <iterator>
#include <cstdint>
#include "fmt/format.h"
//
//template<class ...Args>
//void print(std::string_view formatString, Args&&... args) {
//	//std::cout << std::vformat(formatString, std::make_format_args(args...));
//	//std::vformat_to(std::ostream_iterator<char>(std::cout), formatString, std::make_format_args(args...));
//}
//
//template<class ...Args>
//void println(std::string_view formatString, Args&&... args) {
//	print(formatString, std::forward<Args>(args)...);
//	fmt::format_to(std::ostream_iterator<char>(std::cout), "\n");
//}

template<std::floating_point T, class CharT>
struct fmt::formatter<std::complex<T>, CharT> : fmt::formatter<T, CharT> {

	template<class FormatContext>
	auto format(const std::complex<T>& c, FormatContext& fc) const {
		std::string temp;
		if (!c.imag()) {
			return fmt::formatter<T, CharT>::format(c.real(), fc);
		}
		if (c.real()) {
			fmt::formatter<T, CharT>::format(c.real(), fc);
		}
		if (c.imag()) {
			if (c.imag() > 0.) {
				if (c.real()) fmt::format_to(fc.out(), "+");
				fmt::formatter<T, CharT>::format(c.imag(), fc);
				return fmt::format_to(fc.out(), "i");
			}
			else {
				fmt::format_to(fc.out(), "-");
				fmt::formatter<T, CharT>::format(-c.imag(), fc);
				return fmt::format_to(fc.out(), "i");
			}
		}
		return fmt::format_to(fc.out(), "");
	}
};


namespace Math {
	template <class T, size_t m, size_t n>
	class Matrix;
}
template <class T, size_t m, size_t n>
struct fmt::formatter<Math::Matrix<T, m, n>> : nested_formatter<T> {
	using Parent = nested_formatter<T>;

	auto format(const Math::Matrix<T, m, n>& matrix, format_context& ctx) const {
		auto buffer = std::string{};
		std::vector<size_t> colWidths(matrix.cols());
		for (size_t j = 0; j < matrix.cols(); ++j) {
			colWidths[j] = std::accumulate(matrix.col_begin(j), matrix.col_end(j), size_t{ 0 },
				[&](const auto& max, const auto& el) {
					buffer.clear();
					fmt::format_to(std::back_inserter(buffer), "{}", this->nested(el));
					return std::max(max, buffer.size());
				});
		}
		for (size_t i = 0; i < matrix.rows(); ++i) {
			fmt::format_to(ctx.out(), "| ");
			for (size_t j = 0; j < matrix.cols(); ++j) {
				buffer.clear();
				fmt::format_to(std::back_inserter(buffer), "{}", this->nested(matrix(i, j)));
				fmt::format_to(ctx.out(), "{:{}} ", buffer, colWidths[j]);
			}
			fmt::format_to(ctx.out(), "|\n");
		}
		return ctx.out();
	}
};


template<class T, size_t n, class CharT>
struct fmt::formatter<std::array<T, n>, CharT> : fmt::formatter<T, CharT> {

	template<class FormatContext>
	auto format(const std::array<T, n>& m, FormatContext& fc) const {
		fmt::format_to(fc.out(), "{{");
		if constexpr (n > 0) {
			fmt::formatter<T, CharT>::format(m[0], fc);
			for (size_t i = 1; i < n; ++i) {
				fmt::format_to(fc.out(), ",");
				fmt::formatter<T, CharT>::format(m[i], fc);
			}
		}
		return fmt::format_to(fc.out(), "}}");
	}
};



template<class T, class CharT>
struct fmt::formatter<std::vector<T>, CharT> : fmt::formatter<T, CharT> {

	template<class FormatContext>
	auto format(const std::vector<T>& m, FormatContext& fc) const {
		fmt::format_to(fc.out(), "{{");
		if (m.size() > 0) {
			fmt::formatter<T, CharT>::format(m[0], fc);
			for (size_t i = 1; i < m.size(); ++i) {
				fmt::format_to(fc.out(), ",");
				fmt::formatter<T, CharT>::format(m[i], fc);
			}
		}
		return fmt::format_to(fc.out(), "}}");
	}
};



template<class CharT>
struct fmt::formatter<Q::Binary, CharT> : fmt::formatter<int, CharT> {
	template<class FormatContext>
	auto format(Q::Binary b, FormatContext& fc) const {
		return fmt::format_to(fc.out(), "{}", b.toInt());
	}
};



namespace Q::efficient {
	template<int n>
	class BinaryVector;

	template<int m, int n>
	class BinaryRowMatrix;

	template<int m, int n>
	class BinaryColMatrix;
}

namespace Q {
	template<int n>
	class BinaryPauliOperator;
	
	class Pauli;
}

template<int n, class CharT>
struct fmt::formatter<Q::efficient::BinaryVector<n>, CharT> : fmt::formatter<uint64_t, CharT> {
	template<class FormatContext>
	auto format(Q::efficient::BinaryVector<n> v, FormatContext& fc) const {
		char str[n + 1];
		for (size_t i = 0; i < n; ++i) {
			str[i] = v() & (1ULL << i) ? '1' : '0';
		}
		str[n] = '\0';
		return fmt::format_to(fc.out(), "{}", str);
	}
};


template<int m, int n, class CharT>
struct fmt::formatter<Q::efficient::BinaryRowMatrix<m, n>, CharT> : fmt::formatter<Q::efficient::BinaryVector<n>, CharT> {
	template<class FormatContext>
	auto format(const Q::efficient::BinaryRowMatrix<m, n>& mat, FormatContext& fc) const {
		for (size_t i = 0; i < m; ++i) {
			fmt::format_to(fc.out(), "{}\n", mat[i]);
		}
		return fmt::format_to(fc.out(), "");
	}
};


template<int m, int n, class CharT>
struct fmt::formatter<Q::efficient::BinaryColMatrix<m, n>, CharT> : fmt::formatter<Q::efficient::BinaryVector<m>, CharT> {
	template<class FormatContext>
	auto format(const Q::efficient::BinaryColMatrix<m, n>& mat, FormatContext& fc) const {
		for (size_t i = 0; i < n; ++i) {
			fmt::format_to(fc.out(), "{}\n", mat[i]);
		}
		return fmt::format_to(fc.out(), "");
	}
};

template<int n, class CharT>
struct fmt::formatter<Q::BinaryPauliOperator<n>, CharT> : fmt::formatter<std::string_view, CharT> {
	template<class FormatContext>
	auto format(const Q::BinaryPauliOperator<n>& op, FormatContext& fc) const {
		const auto phase = op.getPhase();
		if (phase != Q::BinaryPhase{ 0 })
			fmt::format_to(fc.out(), "{}", phase.toString());
		return fmt::format_to(fc.out(), "{}", op.toString());
	}
};



template<class T, class CharT>
struct fmt::formatter<std::pair<T,T>, CharT> : fmt::formatter<T, CharT> {
	template<class FormatContext>
	auto format(const std::pair<T,T>& value, FormatContext& fc) const {
		fmt::format_to(fc.out(), "(");
		fmt::formatter<T, CharT>::format(value.first, fc);
		fmt::format_to(fc.out(), ",");
		fmt::formatter<T, CharT>::format(value.second, fc);
		return fmt::format_to(fc.out(), ")");
	}
};

//template<class CharT>
//struct fmt::formatter<Q::BinaryPhase, CharT> : fmt::formatter<std::string_view, CharT> {
//	template<class FormatContext>
//	auto format(const Q::BinaryPhase& phase, FormatContext& fc) const {
//		return fmt::format_to(fc.out(), "{}", phase.toString());
//	}
//};
