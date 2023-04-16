#include "catch2/catch_test_macros.hpp"
#include "catch2/catch_approx.hpp"

#include "pauli.h"


using namespace Q;

TEST_CASE("commutesLocally") {
	REQUIRE(commutesLocally(Pauli{ "II" }, Pauli{ "II" }, 0b11) == true);
	REQUIRE(commutesLocally(Pauli{ "II" }, Pauli{ "II" }, 0b10) == true);
	REQUIRE(commutesLocally(Pauli{ "II" }, Pauli{ "II" }, 0b01) == true);

	REQUIRE(commutesLocally(Pauli{ "XX" }, Pauli{ "YZ" }, 0b11) == true);
	REQUIRE(commutesLocally(Pauli{ "XX" }, Pauli{ "YZ" }, 0b10) == false);
	REQUIRE(commutesLocally(Pauli{ "XX" }, Pauli{ "YZ" }, 0b01) == false);

	REQUIRE(commutesLocally(Pauli{ "XZXXIIX" }, Pauli{ "YIZZXYZ" }, 0b1000111) == false);
}