#include "catch2/catch_test_macros.hpp"
#include "catch2/catch_approx.hpp"

#include "efficient_binary_math.h"


using namespace Q;





TEST_CASE("Efficient Matrix") {
	Math::Matrix<Binary, 4, 4> mat{
		1,0,1,1,
		0,1,0,0,
		0,1,1,1,
		1,0,0,1
	};
	REQUIRE(mat == efficient::BinaryColMatrix<4, 4>{mat}.toMatrix());
	REQUIRE(mat == efficient::BinaryRowMatrix<4, 4>{mat}.toMatrix());
}


TEST_CASE("Efficient binary vector dot product") {
	using Vec = efficient::BinaryVector<8>;
	Vec v{ 0b11001100 };
	Vec w{ 0100000100 };
	REQUIRE(v.dot({ 0b10000100 }) == 0);
	REQUIRE(v.dot({ 0b10000000 }) == 1);
	REQUIRE(v.dot({ 0b00000100 }) == 1);
	REQUIRE(v.dot({ 0b11111111 }) == 0);
	REQUIRE(v.dot({ 0b11111011 }) == 1);
	REQUIRE(v.dot({ 0b00111011 }) == 1);
	REQUIRE(v.dot({ 0b00110011 }) == 0);
}



TEST_CASE("Efficient matrix") {
	Math::Matrix<Binary, 3, 3> m{ 1,0,0,1,1,0,0,1,0 };
	Math::Vector<Binary, 3> v{ 0,0,0 };

	SECTION("Row matrix") {

		efficient::BinaryRowMatrix<3, 3> em{ m };
		REQUIRE(m == em.toMatrix());

		auto vec = efficient::toBitstringInteger<3>(v);

		for (int i = 0; i < 10000; ++i) {
			efficient::BinaryRowMatrix<4, 3> em;
			em.rows[0] = rand();
			em.rows[1] = rand();
			em.rows[2] = rand();
			em.rows[3] = rand();
			efficient::BinaryVector<3> v(rand());
			REQUIRE((em * v).toVector() == em.toMatrix() * v.toVector());
		}
		REQUIRE(efficient::createBinaryVector<3>(em * vec) == m * v);

	}
	SECTION("Col matrix") {

		efficient::BinaryColMatrix<3, 3> em{ m };
		REQUIRE(m == em.toMatrix());


		for (int i = 0; i < 10000; ++i) {
			efficient::BinaryColMatrix<4, 3> em;
			em.cols[0] = rand();
			em.cols[1] = rand();
			em.cols[2] = rand();
			efficient::BinaryVector<4> v(rand());
			REQUIRE((v * em).toVector().transpose() == v.toVector().transpose() * em.toMatrix());
		}

	}
}
