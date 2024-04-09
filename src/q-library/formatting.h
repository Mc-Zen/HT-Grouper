
#pragma once
#include "binary.h"
#include "binary_phase.h"
#include <type_traits>
#include <format>
#include <complex>
#include <iostream>
#include <numeric>
#include <vector>
#include <iterator>
#include <cstdint>

template<class ...Args>
void print(std::string_view formatString, Args&&... args) {
	//std::cout << std::vformat(formatString, std::make_format_args(args...));
	std::vformat_to(std::ostream_iterator<char>(std::cout), formatString, std::make_format_args(args...));
}

template<class ...Args>
void println(std::string_view formatString, Args&&... args) {
	print(formatString, std::forward<Args>(args)...);
	std::format_to(std::ostream_iterator<char>(std::cout), "\n");
}

template<std::floating_point T, class CharT>
struct std::formatter<std::complex<T>, CharT> : std::formatter<T, CharT> {

	template<class FormatContext>
	auto format(const std::complex<T>& c, FormatContext& fc) const {
		std::string temp;
		if (!c.imag()) {
			return std::formatter<T, CharT>::format(c.real(), fc);
		}
		if (c.real()) {
			std::formatter<T, CharT>::format(c.real(), fc);
		}
		if (c.imag()) {
			if (c.imag() > 0.) {
				if (c.real()) std::format_to(fc.out(), "+");
				std::formatter<T, CharT>::format(c.imag(), fc);
				return std::format_to(fc.out(), "i");
			}
			else {
				std::format_to(fc.out(), "-");
				std::formatter<T, CharT>::format(-c.imag(), fc);
				return std::format_to(fc.out(), "i");
			}
		}
		return std::format_to(fc.out(), "");
	}
};


namespace Math {
	template<class T, size_t m, size_t n>
	class Matrix;
}

template<class T, size_t m, size_t n, class CharT>
struct std::formatter<Math::Matrix<T, m, n>, CharT> : std::formatter<T, CharT> {

	std::string formatString;

	// Hack, copy the original format text into a std::string
	constexpr auto parse(format_parse_context& ctx)
	{
		formatString = "{:";
		for (auto iter = std::begin(ctx); iter != std::end(ctx); ++iter) {
			char c = *iter; {
				formatString += c;
			}
			if (c == '}') {
				return iter;
			}
		}
		formatString += '}';
		return std::formatter<T, CharT>::parse(ctx);
	}

	template<class FormatContext>
	auto format(const Math::Matrix<T, m, n>& mat, FormatContext& fc) const {
		//std::array<size_t, n> colWidths;
		//for (size_t j = 0; j < m.cols(); ++j) {
		//	colWidths[j] = std::accumulate(m.col_begin(j), m.col_end(j), 0ULL,
		//		[](const auto& max, const auto& el) { return std::max(max, std::formatted_size("{}", el)); });
		//}
		//
		//for (size_t i = 0; i < m.rows(); ++i) {
		//	std::format_to(fc.out(), "| ");
		//	for (size_t j = 0; j < m.cols(); ++j) {
		//		std::format_to(fc.out(), "{:{}} ", m(i, j), colWidths[j]);
		//	}
		//	std::format_to(fc.out(), "|\n");
		//}
		std::string s;
		std::vector<size_t> colWidths(mat.cols());
		for (size_t j = 0; j < mat.cols(); ++j) {
			colWidths[j] = std::accumulate(mat.col_begin(j), mat.col_end(j), size_t{ 0 },
				[&](const auto& max, const auto& el) {
					s.clear();
					std::vformat_to(std::back_inserter(s), formatString, std::make_format_args(el));
					return std::max(max, s.size());
				});
		}

		for (size_t i = 0; i < mat.rows(); ++i) {
			std::format_to(fc.out(), "| ");
			for (size_t j = 0; j < mat.cols(); ++j) {
				s.clear();
				std::vformat_to(std::back_inserter(s), formatString, std::make_format_args(mat(i, j)));
				std::format_to(fc.out(), "{:{}} ", s, colWidths[j]);
			}
			std::format_to(fc.out(), "|\n");
		}/*
		std::vector<std::string> entries;

		for (size_t i = 0; i < m.rows(); ++i) {
			std::formatter<CharT>::format(m[i], fc);
			std::format_to(fc.out(), "| ");
			for (size_t j = 0; j < m.cols(); ++j) {
				std::format_to(fc.out(), "{:{}} ", m(i, j), colWidths[j]);
				std::format_to(fc.out(), "{:{}} ", m(i, j), colWidths[j]);
			}
			std::format_to(fc.out(), "|\n");
		}*/
		return std::format_to(fc.out(), "");
	}
};


template<class T, size_t n, class CharT>
struct std::formatter<std::array<T, n>, CharT> : std::formatter<T, CharT> {

	template<class FormatContext>
	auto format(const std::array<T, n>& m, FormatContext& fc) const {
		std::format_to(fc.out(), "{{");
		if constexpr (n > 0) {
			std::formatter<T, CharT>::format(m[0], fc);
			for (size_t i = 1; i < n; ++i) {
				std::format_to(fc.out(), ",");
				std::formatter<T, CharT>::format(m[i], fc);
			}
		}
		return std::format_to(fc.out(), "}}");
	}
};



template<class T, class CharT>
struct std::formatter<std::vector<T>, CharT> : std::formatter<T, CharT> {

	template<class FormatContext>
	auto format(const std::vector<T>& m, FormatContext& fc) const {
		std::format_to(fc.out(), "{{");
		if (m.size() > 0) {
			std::formatter<T, CharT>::format(m[0], fc);
			for (size_t i = 1; i < m.size(); ++i) {
				std::format_to(fc.out(), ",");
				std::formatter<T, CharT>::format(m[i], fc);
			}
		}
		return std::format_to(fc.out(), "}}");
	}
};



template<class CharT>
struct std::formatter<Q::Binary, CharT> : std::formatter<int, CharT> {
	template<class FormatContext>
	auto format(Q::Binary b, FormatContext& fc) const {
		return std::format_to(fc.out(), "{}", b.toInt());
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

	class BinaryPhase;
}

template<int n, class CharT>
struct std::formatter<Q::efficient::BinaryVector<n>, CharT> : std::formatter<uint64_t, CharT> {
	template<class FormatContext>
	auto format(Q::efficient::BinaryVector<n> v, FormatContext& fc) const {
		char str[n + 1];
		for (size_t i = 0; i < n; ++i) {
			str[i] = v() & (1ULL << i) ? '1' : '0';
		}
		str[n] = '\0';
		return std::format_to(fc.out(), "{}", str);
	}
};


template<int m, int n, class CharT>
struct std::formatter<Q::efficient::BinaryRowMatrix<m, n>, CharT> : std::formatter<Q::efficient::BinaryVector<n>, CharT> {
	template<class FormatContext>
	auto format(const Q::efficient::BinaryRowMatrix<m, n>& mat, FormatContext& fc) const {
		for (size_t i = 0; i < m; ++i) {
			std::format_to(fc.out(), "{}\n", mat[i]);
		}
		return std::format_to(fc.out(), "");
	}
};


template<int m, int n, class CharT>
struct std::formatter<Q::efficient::BinaryColMatrix<m, n>, CharT> : std::formatter<Q::efficient::BinaryVector<m>, CharT> {
	template<class FormatContext>
	auto format(const Q::efficient::BinaryColMatrix<m, n>& mat, FormatContext& fc) const {
		for (size_t i = 0; i < n; ++i) {
			std::format_to(fc.out(), "{}\n", mat[i]);
		}
		return std::format_to(fc.out(), "");
	}
};

template<int n, class CharT>
struct std::formatter<Q::BinaryPauliOperator<n>, CharT> : std::formatter<std::string_view, CharT> {
	template<class FormatContext>
	auto format(const Q::BinaryPauliOperator<n>& op, FormatContext& fc) const {
		const auto phase = op.getPhase();
		if (phase != Q::BinaryPhase{ 0 })
			std::format_to(fc.out(), "{}", phase.toString());
		return std::format_to(fc.out(), "{}", op.toString());
	}
};



template<class T, class CharT>
struct std::formatter<std::pair<T,T>, CharT> : std::formatter<T, CharT> {
	template<class FormatContext>
	auto format(const std::pair<T,T>& value, FormatContext& fc) const {
		std::format_to(fc.out(), "(");
		std::formatter<T, CharT>::format(value.first, fc);
		std::format_to(fc.out(), ",");
		std::formatter<T, CharT>::format(value.second, fc);
		return std::format_to(fc.out(), ")");
	}
};

//template<class CharT>
//struct std::formatter<Q::BinaryPhase, CharT> : std::formatter<std::string_view, CharT> {
//	template<class FormatContext>
//	auto format(const Q::BinaryPhase& phase, FormatContext& fc) const {
//		return std::format_to(fc.out(), "{}", phase.toString());
//	}
//};
