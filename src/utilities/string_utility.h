#pragma once
#include <string_view>
#include <string>
#include <vector>


namespace Q {
	std::string trim(std::string_view str, std::string_view delims = " ");
	std::string trimFront(std::string_view str, std::string_view delims);
	std::string trimBack(std::string_view str, std::string_view delims);

	std::vector<std::string> split(const std::string_view& s, char delim);

	std::vector<std::string> splitOnce(const std::string& s, char delim);

	std::vector<std::string> split(const std::string_view& str, const std::string_view& delims = " ");


	std::string toLower(std::string_view input);

}