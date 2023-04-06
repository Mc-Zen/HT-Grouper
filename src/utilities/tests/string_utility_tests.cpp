#include "catch2/catch_test_macros.hpp"
#include "catch2/catch_approx.hpp"

#include "string_utility.h"


using namespace Q;


TEST_CASE("String trim") {
	REQUIRE(trim("asd") == "asd");
	REQUIRE(trim(" asd") == "asd");
	REQUIRE(trim("  asd   ") == "asd");
	REQUIRE(trim("  ") == "");
	REQUIRE(trim(" 2") == "2");
	REQUIRE(trim("2 ") == "2");
}

TEST_CASE("String split") {
	//REQUIRE(split("asd") == "asd");

	auto s = split("asd,.dfg,.,.ret,ert,.", ",.");
	REQUIRE(s.size() == 3);
	REQUIRE(s[0] == "asd");
	REQUIRE(s[1] == "dfg");
	REQUIRE(s[2] == "ret,ert");

	s = split("asd,.dfg,.,.ret,ert,.356", ",.");
	REQUIRE(s.size() == 4);
	REQUIRE(s[0] == "asd");
	REQUIRE(s[1] == "dfg");
	REQUIRE(s[2] == "ret,ert");
	REQUIRE(s[3] == "356");
}