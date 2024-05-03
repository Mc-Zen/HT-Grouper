#pragma once

#include <vector>
#include <string>
#include <map>

namespace Q {

	struct CLIArguments {
		std::string path;
		std::vector<std::string> positional;
		std::map<std::string, std::string, std::less<>> options;
	};

	class CLIParseError : public std::runtime_error {
		using std::runtime_error::runtime_error;
	};


	CLIArguments parse_cli_arguments(int argc, char** argv, bool forbidRepeatedArguments = true) {
		CLIArguments args;
		if (argc == 0) return args;

		auto checkRepeatedArgument = [&](bool insertionSuccess, const auto& argument) {
			if (!insertionSuccess && forbidRepeatedArguments) {
				throw CLIParseError(fmt::format("error: the argument '{}' cannnot be used multiple times", argument));
			}
		};

		args.path = argv[0];

		for (int i = 1; i < argc; ++i) {
			std::string argument{ argv[i] };
			if (argument.starts_with("-")) {
				if (auto pos = argument.find("="); pos != std::string::npos) {
					auto [_, success] = args.options.try_emplace(argument.substr(1, pos - 1), argument.substr(pos + 1));
					checkRepeatedArgument(success, argument);
				}
				else {
					auto [_, success] = args.options.try_emplace(argument.substr(1), "");
					checkRepeatedArgument(success, argument);
				}
			}
			else {
				args.positional.emplace_back(argv[i]);
			}
		}
		fmt::println("{}", args.positional);
		for (const auto& [key, val] : args.options) {
			fmt::println("  {}: {}", key, val);
		}
		return args;
	}



}
