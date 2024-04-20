
#pragma once
#include "binary.h"
#include <type_traits>
#include <complex>
#include <iostream>
#include <numeric>
#include <vector>

namespace std{
	template <class T, class _CharT>
	struct formatter{};

	template<class ...Args>
	auto format(Args...args) {
		return std::string{};
	}
	template<class ...Args>
	auto format_to(Args...args) {
		return std::string{};
	}
}

template<class ...Args>
void print(std::string_view formatString, Args&&... args) {
	//std::cout << std::vformat(formatString, std::make_format_args(args...));
	//std::vformat_to(std::ostream_iterator<char>(std::cout), formatString, std::make_format_args(args...));
}

template<class ...Args>
void println(std::string_view formatString, Args&&... args) {
	print(formatString, std::forward<Args>(args)...);
	//std::format_to(std::ostream_iterator<char>(std::cout), "\n");
}