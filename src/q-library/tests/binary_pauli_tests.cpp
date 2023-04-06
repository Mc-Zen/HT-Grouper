#include "catch2/catch_test_macros.hpp"
#include "catch2/catch_approx.hpp"

#include "binary_pauli.h"
#include "formatting.h"


using namespace Q;



TEST_CASE("Binary Phase") {
	BinaryPhase p{};
	REQUIRE(p.toInt() == 0);
	p += 5;
	REQUIRE(p.toInt() == 1);
	p += 5;
	REQUIRE(p.toInt() == 2);
	p -= 5;
	REQUIRE(p.toInt() == 1);
	p += BinaryPhase{ 2 };
	REQUIRE(p.toInt() == 3);
	p -= BinaryPhase{ 2 };
	REQUIRE(p.toInt() == 1);
	REQUIRE((BinaryPhase{ 2 } + BinaryPhase{ 5 }).toInt() == 3);
	REQUIRE((BinaryPhase{ 2 } - BinaryPhase{ 5 }).toInt() == 1);
}


TEST_CASE("Binary Pauli") {
	BinaryPauliOperator<4> op{ "IXYZ" };

	SECTION("No phase") {
		op = BinaryPauliOperator < 4>{ "IXYZ" };
		REQUIRE(op.getPhase() == BinaryPhase{ 0 });
		REQUIRE(op.getXZPhase() == BinaryPhase{ 1 });
		REQUIRE(std::format("{}", op) == "IXYZ");
	}
	SECTION("Phase i") {
		op = BinaryPauliOperator < 4>{ "iIXYZ" };
		REQUIRE(op.getPhase() == BinaryPhase{ 1 });
		REQUIRE(op.getXZPhase() == BinaryPhase{ 2 });
		REQUIRE(std::format("{}", op) == "iIXYZ");
	}
	SECTION("Phase -i") {
		op = BinaryPauliOperator < 4>{ "-iIXYZ" };
		REQUIRE(op.getPhase() == BinaryPhase{ 3 });
		REQUIRE(op.getXZPhase() == BinaryPhase{ 0 });
		REQUIRE(std::format("{}", op) == "-iIXYZ");
	}
	SECTION("Phase -") {
		op = BinaryPauliOperator < 4>{ "-IXYZ" };
		REQUIRE(op.getPhase() == BinaryPhase{ 2 });
		REQUIRE(op.getXZPhase() == BinaryPhase{ 3 });
		REQUIRE(std::format("{}", op) == "-IXYZ");
	}
	SECTION("Phase +") {
		op = BinaryPauliOperator < 4>{ "+IXYZ" };
		REQUIRE(op.getPhase() == BinaryPhase{ 0 });
		REQUIRE(op.getXZPhase() == BinaryPhase{ 1 });
		REQUIRE(std::format("{}", op) == "IXYZ");
	}

	REQUIRE(op[0] == BinaryPauli::I);
	REQUIRE(op[1] == BinaryPauli::X);
	REQUIRE(op[2] == BinaryPauli::Y);
	REQUIRE(op[3] == BinaryPauli::Z);
}

TEST_CASE("Binary Pauli Single Qubit Clifford Operations") {

	BinaryPauliOperator<1> i{ "I" };
	BinaryPauliOperator<1> x{ "X" };
	BinaryPauliOperator<1> y{ "Y" };
	BinaryPauliOperator<1> z{ "Z" };

	SECTION("X") {
		Clifford::x(i, 0);
		Clifford::x(x, 0);
		Clifford::x(y, 0);
		Clifford::x(z, 0);
		REQUIRE(i[0] == BinaryPauli::I);
		REQUIRE(i.getPhase() == 0);
		REQUIRE(x[0] == BinaryPauli::X);
		REQUIRE(x.getPhase() == 0);
		REQUIRE(y[0] == BinaryPauli::Y);
		REQUIRE(y.getPhase() == 2);
		REQUIRE(z[0] == BinaryPauli::Z);
		REQUIRE(z.getPhase() == 2);
	}

	SECTION("Y") {
		Clifford::y(i, 0);
		Clifford::y(x, 0);
		Clifford::y(y, 0);
		Clifford::y(z, 0);
		REQUIRE(i[0] == BinaryPauli::I);
		REQUIRE(i.getPhase() == 0);
		REQUIRE(x[0] == BinaryPauli::X);
		REQUIRE(x.getPhase() == 2);
		REQUIRE(y[0] == BinaryPauli::Y);
		REQUIRE(y.getPhase() == 0);
		REQUIRE(z[0] == BinaryPauli::Z);
		REQUIRE(z.getPhase() == 2);
	}

	SECTION("Z") {
		Clifford::z(i, 0);
		Clifford::z(x, 0);
		Clifford::z(y, 0);
		Clifford::z(z, 0);
		REQUIRE(i[0] == BinaryPauli::I);
		REQUIRE(i.getPhase() == 0);
		REQUIRE(x[0] == BinaryPauli::X);
		REQUIRE(x.getPhase() == 2);
		REQUIRE(y[0] == BinaryPauli::Y);
		REQUIRE(y.getPhase() == 2);
		REQUIRE(z[0] == BinaryPauli::Z);
		REQUIRE(z.getPhase() == 0);
	}

	SECTION("H") {
		Clifford::h(i, 0);
		Clifford::h(x, 0);
		Clifford::h(y, 0);
		Clifford::h(z, 0);
		REQUIRE(i[0] == BinaryPauli::I);
		REQUIRE(i.getPhase() == 0);
		REQUIRE(x[0] == BinaryPauli::Z);
		REQUIRE(x.getPhase() == 0);
		REQUIRE(y[0] == BinaryPauli::Y);
		REQUIRE(y.getPhase() == 2);
		REQUIRE(z[0] == BinaryPauli::X);
		REQUIRE(z.getPhase() == 0);
	}

	SECTION("S") {
		Clifford::s(i, 0);
		Clifford::s(x, 0);
		Clifford::s(y, 0);
		Clifford::s(z, 0);
		REQUIRE(i[0] == BinaryPauli::I);
		REQUIRE(i.getPhase() == 0);
		REQUIRE(x[0] == BinaryPauli::Y);
		REQUIRE(x.getPhase() == 0);
		REQUIRE(y[0] == BinaryPauli::X);
		REQUIRE(y.getPhase() == 2);
		REQUIRE(z[0] == BinaryPauli::Z);
		REQUIRE(z.getPhase() == 0);
	}

	SECTION("SDG") {
		Clifford::sdg(i, 0);
		Clifford::sdg(x, 0);
		Clifford::sdg(y, 0);
		Clifford::sdg(z, 0);
		REQUIRE(i[0] == BinaryPauli::I);
		REQUIRE(i.getPhase() == 0);
		REQUIRE(x[0] == BinaryPauli::Y);
		REQUIRE(x.getPhase() == 2);
		REQUIRE(y[0] == BinaryPauli::X);
		REQUIRE(y.getPhase() == 0);
		REQUIRE(z[0] == BinaryPauli::Z);
		REQUIRE(z.getPhase() == 0);
	}

	SECTION("HS") {
		Clifford::hs(i, 0);
		Clifford::hs(x, 0);
		Clifford::hs(y, 0);
		Clifford::hs(z, 0);
		REQUIRE(i[0] == BinaryPauli::I);
		REQUIRE(i.getPhase() == 0);
		REQUIRE(x[0] == BinaryPauli::Y);
		REQUIRE(x.getPhase() == 2);
		REQUIRE(y[0] == BinaryPauli::Z);
		REQUIRE(y.getPhase() == 2);
		REQUIRE(z[0] == BinaryPauli::X);
		REQUIRE(z.getPhase() == 0);
	}

	SECTION("SH") {
		Clifford::sh(i, 0);
		Clifford::sh(x, 0);
		Clifford::sh(y, 0);
		Clifford::sh(z, 0);
		REQUIRE(i[0] == BinaryPauli::I);
		REQUIRE(i.getPhase() == 0);
		REQUIRE(x[0] == BinaryPauli::Z);
		REQUIRE(x.getPhase() == 0);
		REQUIRE(y[0] == BinaryPauli::X);
		REQUIRE(y.getPhase() == 0);
		REQUIRE(z[0] == BinaryPauli::Y);
		REQUIRE(z.getPhase() == 0);
	}

	SECTION("HSH") {
		Clifford::hsh(i, 0);
		Clifford::hsh(x, 0);
		Clifford::hsh(y, 0);
		Clifford::hsh(z, 0);
		REQUIRE(i[0] == BinaryPauli::I);
		REQUIRE(i.getPhase() == 0);
		REQUIRE(x[0] == BinaryPauli::X);
		REQUIRE(x.getPhase() == 0);
		REQUIRE(y[0] == BinaryPauli::Z);
		REQUIRE(y.getPhase() == 0);
		REQUIRE(z[0] == BinaryPauli::Y);
		REQUIRE(z.getPhase() == 2);
	}
}

TEST_CASE("Binary Pauli 2-Qubit Clifford Operations") {

	BinaryPauliOperator<2> II{ "II" };
	BinaryPauliOperator<2> IX{ "IX" };
	BinaryPauliOperator<2> IY{ "IY" };
	BinaryPauliOperator<2> IZ{ "IZ" };
	BinaryPauliOperator<2> XI{ "XI" };
	BinaryPauliOperator<2> XX{ "XX" };
	BinaryPauliOperator<2> XY{ "XY" };
	BinaryPauliOperator<2> XZ{ "XZ" };
	BinaryPauliOperator<2> YI{ "YI" };
	BinaryPauliOperator<2> YX{ "YX" };
	BinaryPauliOperator<2> YY{ "YY" };
	BinaryPauliOperator<2> YZ{ "YZ" };
	BinaryPauliOperator<2> ZI{ "ZI" };
	BinaryPauliOperator<2> ZX{ "ZX" };
	BinaryPauliOperator<2> ZY{ "ZY" };
	BinaryPauliOperator<2> ZZ{ "ZZ" };

	SECTION("CX") {
		Clifford::cx(II, 0, 1);
		Clifford::cx(IX, 0, 1);
		Clifford::cx(IY, 0, 1);
		Clifford::cx(IZ, 0, 1);
		Clifford::cx(XI, 0, 1);
		Clifford::cx(XX, 0, 1);
		Clifford::cx(XY, 0, 1);
		Clifford::cx(XZ, 0, 1);
		Clifford::cx(YI, 0, 1);
		Clifford::cx(YX, 0, 1);
		Clifford::cx(YY, 0, 1);
		Clifford::cx(YZ, 0, 1);
		Clifford::cx(ZI, 0, 1);
		Clifford::cx(ZX, 0, 1);
		Clifford::cx(ZY, 0, 1);
		Clifford::cx(ZZ, 0, 1);
		REQUIRE(II == "II");
		REQUIRE(IX == "IX");
		REQUIRE(IY == "ZY");
		REQUIRE(IZ == "ZZ");

		REQUIRE(XI == "XX");
		REQUIRE(XX == "XI");
		REQUIRE(XY == "YZ");
		REQUIRE(XZ == "-YY");

		REQUIRE(YI == "YX");
		REQUIRE(YX == "YI");
		REQUIRE(YY == "-XZ");
		REQUIRE(YZ == "XY");

		REQUIRE(ZI == "ZI");
		REQUIRE(ZX == "ZX");
		REQUIRE(ZY == "IY");
		REQUIRE(ZZ == "IZ");
	}

	SECTION("CZ") {
		Clifford::cz(II, 0, 1);
		Clifford::cz(IX, 0, 1);
		Clifford::cz(IY, 0, 1);
		Clifford::cz(IZ, 0, 1);
		Clifford::cz(XI, 0, 1);
		Clifford::cz(XX, 0, 1);
		Clifford::cz(XY, 0, 1);
		Clifford::cz(XZ, 0, 1);
		Clifford::cz(YI, 0, 1);
		Clifford::cz(YX, 0, 1);
		Clifford::cz(YY, 0, 1);
		Clifford::cz(YZ, 0, 1);
		Clifford::cz(ZI, 0, 1);
		Clifford::cz(ZX, 0, 1);
		Clifford::cz(ZY, 0, 1);
		Clifford::cz(ZZ, 0, 1);
		REQUIRE(II == "II");
		REQUIRE(IX == "ZX");
		REQUIRE(IY == "ZY");
		REQUIRE(IZ == "IZ");

		REQUIRE(XI == "XZ");
		REQUIRE(XX == "YY");
		REQUIRE(XY == "-YX");
		REQUIRE(XZ == "XI");

		REQUIRE(YI == "YZ");
		REQUIRE(YX == "-XY");
		REQUIRE(YY == "XX");
		REQUIRE(YZ == "YI");

		REQUIRE(ZI == "ZI");
		REQUIRE(ZX == "IX");
		REQUIRE(ZY == "IY");
		REQUIRE(ZZ == "ZZ");
	}

}
