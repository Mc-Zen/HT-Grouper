#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>


#define MATRIX_EXCEPTIONS
#include "matrix.h"

#include <iostream>


using Catch::Approx;
using namespace Math;


template<std::random_access_iterator It>
void checkRandomAccessIterator(It) {
}


template<std::bidirectional_iterator It>
void checkBirectionalIterator(It) {
}

template<class T, Index m, Index n>
void testSizeImpl() {
	Matrix<float, m, n> mat;
	REQUIRE(m == mat.rows());
	REQUIRE(n == mat.cols());
	REQUIRE(n * m == mat.size());
}


TEST_CASE("Dynamic Constructor with deduction guide") {
	Matrix mat{ 1,2,{4,5,6} };
	static_assert(std::same_as<decltype(mat)::value_type, int>);
	REQUIRE(mat(0, 0) == 4);
}
TEST_CASE("Dynamic Constructor") {
	SECTION("Dimensions") {
		Matrix<int> mat{ Shape{3,2} };
		REQUIRE(mat.rows() == 3);
		REQUIRE(mat.cols() == 2);
	}

	SECTION("m,n") {
		Matrix<int> mat{ 3,2 };
		REQUIRE(mat.rows() == 3);
		REQUIRE(mat.cols() == 2);
	}


	SECTION("Value") {
		Matrix<int> mat{ 3,2,4 };
		REQUIRE(mat.rows() == 3);
		REQUIRE(mat.cols() == 2);
		REQUIRE(mat(0, 0) == 4);
	}


	SECTION("Initializer list") {
		Matrix<int> mat{ 3,2, {3,4,5,6} };
		REQUIRE(mat.rows() == 3);
		REQUIRE(mat.cols() == 2);
		REQUIRE(mat(0, 0) == 3);
		REQUIRE(mat(0, 1) == 4);
		REQUIRE(mat(1, 0) == 5);
		REQUIRE(mat(1, 1) == 6);
		REQUIRE(mat(2, 0) == 0);
		REQUIRE(mat(2, 1) == 0);
	}


	SECTION("vector copy") {
		std::vector<int> data{ 3, 4, 5, 6 };
		Matrix<int> mat{ 3,2, data };
		REQUIRE(mat.rows() == 3);
		REQUIRE(mat.cols() == 2);
		REQUIRE(mat(0, 0) == 3);
		REQUIRE(mat(0, 1) == 4);
		REQUIRE(mat(1, 0) == 5);
		REQUIRE(mat(1, 1) == 6);
		REQUIRE(mat(2, 0) == 0);
		REQUIRE(mat(2, 1) == 0);
	}

	SECTION("vector move") {
		std::vector<int> data{ 3, 4, 5, 6 };
		Matrix<int> mat{ 3,2,  std::move(data) };
		REQUIRE(mat.rows() == 3);
		REQUIRE(mat.cols() == 2);
		REQUIRE(mat(0, 0) == 3);
		REQUIRE(mat(0, 1) == 4);
		REQUIRE(mat(1, 0) == 5);
		REQUIRE(mat(1, 1) == 6);
		REQUIRE(mat(2, 0) == 0);
		REQUIRE(mat(2, 1) == 0);
	}


	SECTION("From Matrix view") {
		Matrix<int> mat{ 4,5,2 };
		Matrix mat2 = mat.row(2);
		REQUIRE(mat2.rows() == 1);
		REQUIRE(mat2.cols() == 5);
		REQUIRE(mat(0, 0) == 2);
		REQUIRE(mat(0, 1) == 2);
		REQUIRE(mat(0, 2) == 2);
		REQUIRE(mat(0, 3) == 2);
		REQUIRE(mat(0, 4) == 2);
	}

}


TEST_CASE("Dynamic arithmetic") {
	Matrix<int> a{ 2,3 };
	Matrix<int> b{ 2,3 };
	a + b;
	b.resize(2, 4);

	bool exceptionHappened{ false };
	try {
		a + b;
	}
	catch (Matrix_shape_error&) {
		exceptionHappened = true;
	}
	REQUIRE(exceptionHappened);

	b.resize(3, 5);
	a = a * b;
	exceptionHappened = false;
	try {
		a* b;
	}
	catch (Matrix_shape_error&) {
		exceptionHappened = true;
	}
	REQUIRE(exceptionHappened);
}


TEST_CASE("Dynamic transpose") {
	Matrix<int> a{ 2,3 };
	a = a.transpose();
	REQUIRE(a.rows() == 3);
	REQUIRE(a.cols() == 2);
}

TEST_CASE("Dynamic dot product") {
	Matrix<int> a{ 2,1, {3,4} };
	Matrix<int> b{ 2,1, {6,7} };
	REQUIRE(a.dot(b) == 18 + 28);
	b.resize(2, 2);
	bool exceptionHappened{ false };
	try {
		a.dot(b);
	}
	catch (Matrix_shape_error&) {
		exceptionHappened = true;
	}
	REQUIRE(exceptionHappened);
}


TEST_CASE("Dynamic diag") {
	SECTION("Vector") {
		Matrix vec(4, 1, { 3,4,2,1 });
		auto mat = diag(vec);
		REQUIRE(mat.rows() == 4);
		REQUIRE(mat.cols() == 4);
		REQUIRE(mat(0, 0) == 3);
		REQUIRE(mat(1, 1) == 4);
		REQUIRE(mat(2, 2) == 2);
		REQUIRE(mat(3, 3) == 1);
	}

	SECTION("Initializer list") {
		auto mat = diag({ 3,4,2,1 });
		REQUIRE(mat.rows() == 4);
		REQUIRE(mat.cols() == 4);
		REQUIRE(mat(0, 0) == 3);
		REQUIRE(mat(1, 1) == 4);
		REQUIRE(mat(2, 2) == 2);
		REQUIRE(mat(3, 3) == 1);
	}
}


TEST_CASE("Dynamic antidiag") {
	SECTION("Vector") {
		Matrix vec(4, 1, { 3,4,2,1 });
		auto mat = antidiag(vec);
		REQUIRE(mat.rows() == 4);
		REQUIRE(mat.cols() == 4);
		REQUIRE(mat(3, 0) == 3);
		REQUIRE(mat(2, 1) == 4);
		REQUIRE(mat(1, 2) == 2);
		REQUIRE(mat(0, 3) == 1);
	}

	SECTION("Initializer list") {
		auto mat = antidiag({ 3,4,2,1 });
		REQUIRE(mat.rows() == 4);
		REQUIRE(mat.cols() == 4);
		REQUIRE(mat(3, 0) == 3);
		REQUIRE(mat(2, 1) == 4);
		REQUIRE(mat(1, 2) == 2);
		REQUIRE(mat(0, 3) == 1);
	}
}

TEST_CASE("Random Acess iterator concept") {
	Matrix<float, 4, 3> mat;
	checkRandomAccessIterator(mat.begin());
	checkRandomAccessIterator(mat.col_begin(0));
	checkBirectionalIterator(mat.block(1, 1, 1, 1).begin());
}


TEST_CASE("Size") {
	testSizeImpl<float, 3, 4>();
	testSizeImpl<float, 3, 9>();
	testSizeImpl<float, 3, 21>();
	testSizeImpl<float, 1, 1>();
	testSizeImpl<float, 2, 1>();
	// Matrix<float, 4, 5> a
}

TEST_CASE("EmptyConstructor") {
	Matrix<float, 4, 3> mat;
	for (const auto& el : mat) {
		REQUIRE(el == 0.f);
	}
}

TEST_CASE("CopyConstructor") {
	Matrix<float, 4, 3> mat(3);
	Matrix<float, 4, 3> copy(mat);
	for (const auto& el : copy) {
		REQUIRE(el == 3.f);
	}
}

TEST_CASE("MoveConstructor") {
	Matrix<float, 4, 3> mat(3);
	Matrix<float, 4, 3> moved(std::move(mat));
	for (const auto& el : moved) {
		REQUIRE(el == 3.f);
	}
}

TEST_CASE("CopyAssignment") {
	Matrix<float, 4, 3> mat(3);
	Matrix<float, 4, 3> copy = mat;
	for (const auto& el : copy) {
		REQUIRE(el == 3.f);
	}
}

TEST_CASE("MoveAssignment") {
	Matrix<float, 4, 3> mat(3);
	Matrix<float, 4, 3> moved = std::move(mat);
	for (const auto& el : moved) {
		REQUIRE(el == 3.f);
	}
}


TEST_CASE("OneValueConstructor") {
	Matrix<float, 4, 3> mat(4);
	for (const auto& el : mat) {
		REQUIRE(el == 4.f);
	}
}

TEST_CASE("InitializerListConstructor1") {
	Matrix<float, 4, 3> mat{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
	float count = 0;
	for (const auto& el : mat) {
		REQUIRE(el == count++);
	}
}

TEST_CASE("InitializerListConstructor2") {
	Matrix<float, 4, 3> mat = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
	float count = 0;
	for (const auto& el : mat) {
		REQUIRE(el == count++);
	}
}

TEST_CASE("InitializerListConstructorTooManyValues") {
	Matrix<float, 4, 3> mat = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 };
	float count = 0;
	for (const auto& el : mat) {
		REQUIRE(el == count++);
	}
}

TEST_CASE("InitializerListConstructorTooFewValues") {
	Matrix<float, 4, 3> mat = { 0, 1, 2, 3, 4, 5 };
	float count = 0;
	for (const auto& el : mat) {
		REQUIRE(el == (count < 6 ? count++ : 0.f));
	}
}

TEST_CASE("ArrayCopyConstructor") {
	std::array<float, 12> arr{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
	Matrix<float, 4, 3> mat{ arr };
	float count = 0;
	for (const auto& el : mat) {
		REQUIRE(el == count++);
	}
}

TEST_CASE("ArrayMoveConstructor") {
	std::array<float, 12> arr{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
	Matrix<float, 4, 3> mat{ std::move(arr) };
	float count = 0;
	for (const auto& el : mat) {
		REQUIRE(el == count++);
	}
}

TEST_CASE("MatrixViewConstructor") {
	Matrix<float, 4, 3> mat(4.f);
	Matrix<float, 3, 2> mat2{ mat.block(0, 0, 3, 2) };
	Vector<float, 4> vec{ mat.col(1) };

	for (const auto& el : mat2)
		REQUIRE(el == 4.f);
	for (const auto& el : vec)
		REQUIRE(el == 4.f);
	bool exceptionHappened{ false };
	try {
		Vector<float, 3> vec1{ mat.col(1) };
	}
	catch (Matrix_block_domain_error&) {
		exceptionHappened = true;
	}
	REQUIRE(exceptionHappened);
}

TEST_CASE("IdentityMatrixFactory") {
	constexpr int d = 7;
	auto mat = Matrix<float, d, d>::identity();
	for (Index i = 0; i < d; i++)
		for (Index j = 0; j < d; j++)
			REQUIRE(mat(i, j) == (i == j ? 1.f : 0.f));
}

TEST_CASE("ZeroMatrixFactory") {
	auto mat = Matrix<float, 7, 9>::zero();
	for (const auto& c : mat)
		REQUIRE(c == 0.f);
}

TEST_CASE("DiagMatrixFactory") {
	auto mat = diag<double, 3>({ 9., 10., 11. });
	double num = 9.;
	for (Index i = 0; i < mat.rows(); i++)
		for (Index j = 0; j < mat.rows(); j++)
			REQUIRE(mat(i, j) == (i == j ? num++ : 0.));
}

TEST_CASE("AntiDiagMatrixFactory") {
	auto mat = antidiag({ 9., 10., 11. });
	double num = 9.;
	for (Index i = 0; i < mat.rows(); i++)
		for (Index j = 0; j < mat.rows(); j++)
			REQUIRE(mat(mat.rows() - i - 1, j) == (i == j ? num++ : 0.));
}

TEST_CASE("ParenthesesAccessOperator") {
	std::array<float, 12> arr{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
	Matrix<float, 4, 3> mat{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
	auto it = mat.begin();
	for (Index i = 0; i < 4; i++) {
		for (Index j = 0; j < 3; j++) {
			REQUIRE(*(it++) == mat(i, j));
		}
	}
}

TEST_CASE("AtAccessOperator") {
	std::array<float, 12> arr{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
	Matrix<float, 4, 3> mat{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
	auto it = mat.begin();
	for (Index i = 0; i < 4; i++) {
		for (Index j = 0; j < 3; j++) {
			REQUIRE(*(it++) == mat.at(i, j));
		}
	}
}

TEST_CASE("AtAccessOperatorBoundsException") {
	Matrix<float, 3, 4> mat;
	bool exceptionHappened{ false };
	try {
		mat.at(3, 12);
	}
	catch (std::out_of_range&) {
		exceptionHappened = true;
	}
	REQUIRE(exceptionHappened == true);
}

TEST_CASE("ConstIterator") {
	Matrix<float, 4, 3> mat = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
	for (auto& a : mat) {
		a = -1;
	}
	for (const auto& a : mat) {
		REQUIRE(-1.f == a);
	}
}

TEST_CASE("Iterator") {
	Matrix<float, 4, 3> mat = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
	for (auto& a : mat) {
		a++;
	}
	float count = 1;
	for (const auto& a : mat) {
		REQUIRE(count++ == a);
	}
}

TEST_CASE("ReverseIterator") {
	Matrix<float, 4, 3> mat = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
	for (auto it = mat.rbegin(); it != mat.rend(); ++it) {
		(*it)++;
	}
	float count = 12;
	for (auto it = mat.rbegin(); it != mat.rend(); ++it) {
		REQUIRE(count-- == *it);
	}
}

TEST_CASE("ColIterator") {
	Matrix<float, 3, 4> mat{
		0, 3, 6, 9,
		1, 4, 7, 10,
		2, 5, 8, 11
	};
	float count = 0;
	for (Index col = 0; col < mat.cols(); col++) {
		for (auto it = mat.col_begin(col); it != mat.col_end(col); it++) {
			REQUIRE(count++ == *it);
		}
	}
}

TEST_CASE("RowIterator") {
	Matrix<float, 3, 4> mat{
		0, 1, 2, 3,
		4, 5, 6, 7,
		8, 9, 10, 11
	};
	float count = 0;
	for (Index row = 0; row < mat.rows(); row++) {
		for (auto it = mat.row_begin(row); it != mat.row_end(row); it++) {
			REQUIRE(count++ == *it);
		}
	}
}

TEST_CASE("Arithmetic") {
	Matrix<float, 3, 4> mat1(2);
	Matrix<float, 3, 4> mat2(-7);
	auto mat3 = mat1 + mat2;
	for (const auto& el : mat3)
		REQUIRE(-5.f == el);
	mat3 = mat1 * 2.f - mat2;
	for (const auto& el : mat3)
		REQUIRE(11.f == el);
	mat3 = -mat1;
	for (const auto& el : mat3)
		REQUIRE(-2.f == el);
	Matrix<int, 3, 4> mat4(8);
	auto mat5 = mat4 / 2;
	for (const auto& el : mat5)
		REQUIRE(4 == el);
}

TEST_CASE("Transposition") {
	Matrix<float, 3, 4> mat;
	float index = 0;
	for (auto& c : mat)
		c = index++;
	auto transposed = mat.transpose();
	index = 0;
	for (Index j = 0; j < transposed.cols(); ++j)
		for (Index i = 0; i < transposed.rows(); ++i)
			REQUIRE(index++ == transposed(i, j));
}

TEST_CASE("MatrixMultiplication") {
	Matrix<float, 3, 4> mat1;
	Matrix<float, 4, 8> mat2;
	float index = 0;
	for (auto& c : mat1)
		c = index++;
	for (auto& c : mat2)
		c = index++;
	float numbers[] = { 178, 546, 914 };
	float differences[] = { 6, 22, 38 };
	auto product = mat1 * mat2;
	for (Index i = 0; i < product.rows(); ++i) {
		for (const auto& c : product.row(i))
			REQUIRE((numbers[i] += differences[i]) == c);
	}
}

TEST_CASE("Comparison") {
	Matrix<float, 3, 4> mat1(2);
	Matrix<float, 3, 4> mat2(-7);
	REQUIRE((mat1 == mat1) == true);
	REQUIRE((mat1 != mat1) == false);
	REQUIRE((mat1 == mat2) == false);
	REQUIRE((mat1 != mat2) == true);
	REQUIRE(mat1 == mat1);
	REQUIRE(mat1 == mat2 + 9.f);
}

TEST_CASE("Assignment") {
	Matrix<float, 3, 4> mat1(2);
	Matrix<float, 3, 4> mat2(-7);
	REQUIRE(mat1 != mat2);
	mat1 = mat2;
	REQUIRE(mat1 == mat2);
}

TEST_CASE("VectorAccess") {
	Vector<float, 7> vec(2);
	for (int i = 0; i < 7; i++) {
		REQUIRE(2.f == vec[i]);
		vec[i] = -10.f;
	}
	for (const auto& el : vec) {
		REQUIRE(-10.f == el);
	}
	Vector<float, 1> v1(1);
	v1.x() = 9.f;
	REQUIRE(9.f == v1.x());
	Vector<float, 2> v2(1);
	v2.x() = 9.f;
	v2.y() = 5.f;
	REQUIRE(5.f == v2.y());
	Vector<float, 3> v3(1);
	v3.x() = 9.f;
	v3.y() = 5.f;
	v3.z() = 3.f;
	REQUIRE(3.f == v3.z());

	Vector<float, 4> v4(1);
	v4.x() = 9.f;
	v4.y() = 5.f;
	v4.z() = 3.f;
	v4.w() = -1.f;
	REQUIRE(-1.f == v4.w());
}

TEST_CASE("VectorNorm") {
	Vector<float, 4> vec(3);
	REQUIRE(vec.norm() == Approx(6.f));

	Vector<float, 8> vec2{ 1, 2, 3, 4, 5, 6, 10 };
	REQUIRE(vec2.normalize().norm() == Approx(1.0f));
}

TEST_CASE("VectorInnerProduct") {
	Vector<float, 3> vec1{ 2, 3, -5 };
	Vector<float, 3> vec2{ -34, 2, 99 };
	REQUIRE(-557.f == vec1 * vec2);
}

TEST_CASE("VectorDistance") {
	Vector<float, 3> vec1{ 0, 0, 0 };
	Vector<float, 3> vec2{ 2, 2, 2 };
	REQUIRE(distance(vec1, vec2) == Approx(std::sqrtf(12.f)));
	RowVector<float, 3> rvec1{ 0, 0, 0 };
	RowVector<float, 3> rvec2{ 2, 2, 2 };
	REQUIRE(distance(rvec1, rvec2) == Approx(std::sqrtf(12.f)));
}

TEST_CASE("Matrix1x1CastToT") {
	Matrix<float, 1, 1> mat(1);
	float value = mat;
	REQUIRE(1.f == value);
}

TEST_CASE("MatrixFill") {
	Matrix<float, 4, 2> mat(1);
	mat.fill(2.f);
	for (const auto& c : mat)
		REQUIRE(c == 2.f);
}

TEST_CASE("MatrixSwap") {
	Matrix<float, 4, 2> mat1(1);
	Matrix<float, 4, 2> mat2(3);
	mat1.swap(mat2);
	std::swap(mat1, mat2);
	for (const auto& c : mat1)
		REQUIRE(c == 1.f);
	for (const auto& c : mat2)
		REQUIRE(c == 3.f);
}

TEST_CASE("CastValueTypeToT") {
	Matrix<double, 4, 4> mat(2.2);
	auto mat_int = mat.cast<int>();
	auto mat_bool = mat.cast<bool>();
	for (const auto& c : mat_int)
		REQUIRE(c == 2);
	for (const auto& c : mat_bool)
		REQUIRE(c == true);
}

TEST_CASE("MatrixBlockIterator") {
	Matrix<float, 4, 4> mat{ 11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34, 41, 42, 43, 44 };
	std::array<float, 6> expected{ 22, 23, 32, 33, 42, 43 };
	int index = 0;
	for (auto& c : mat.block(1, 1, 3, 3))
		c++;
	for (const auto& c : mat.block(1, 1, 3, 2)) {
		REQUIRE(expected[index] + 1 == c);
		index++;
	}

	auto block = mat.block(1, 1, 3, 2);
	for (auto it = block.rbegin(); it != block.rend(); it++)
		(*it)--;
	index = 6;
	for (auto it = block.rbegin(); it != block.rend(); it++) {
		--index;
		REQUIRE(expected[index] == *it);
	}
}

TEST_CASE("MatrixRow") {
	Matrix<float, 4, 4> mat{ 11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34, 41, 42, 43, 44 };
	for (size_t i = 0; i < mat.rows(); i++) {
		for (auto& c : mat.row(i))
			c--;
		int index = 0;
		for (const auto& c : mat.row(i))
			REQUIRE((i + 1) * 10.f + (index++) == c);
	}
}

TEST_CASE("MatrixCol") {
	Matrix<float, 4, 4> mat{ 11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34, 41, 42, 43, 44 };
	for (size_t i = 0; i < mat.cols(); i++) {
		for (auto& c : mat.col(i))
			c--;
		int index = 0;
		for (const auto& c : mat.col(i))
			REQUIRE(i + 10.f * (++index) == c);
	}
}

TEST_CASE("MatrixView") {
	Matrix<float, 6, 7> mat1(1);
	Matrix<float, 4, 4> mat2(3);

	mat1.block(1, 2, 2, 3) = mat2.block(1, 1, 2, 3);
	std::cout << mat1 << "\n";
	for (size_t i = 0; i < mat1.rows(); ++i) {
		for (size_t j = 0; j < mat1.cols(); ++j) {
			REQUIRE(mat1(i, j) == ((i > 0 && i < 3 && j > 1 && j < 5) ? 3.f : 1.f));
		}
	}
	float index = 0;
	for (auto& c : mat1)
		c = index++;
	index = 0;
	for (auto& c : mat2)
		c = index--;


	mat1.block(0, 2, 2, 3) = mat2.block(1, 1, 2, 3);
	mat2.row(3) = mat2.row(2);


	std::cout << mat1 << "\n";
	std::cout << mat2 << "\n";
}

TEST_CASE("CopyBlockToMatrix") {
	Matrix<float, 6, 7> mat1(1);
	Matrix<float, 4, 4> mat2(3);
	mat2 = mat1.block(0, 0, 4, 4);
	for (const auto& c : mat2)
		REQUIRE(c == 1.f);
	bool exceptionHappened{ false };
	try {
		mat2 = mat1.block(0, 0, 4, 3);
	}
	catch (std::exception&) {
		exceptionHappened = true;
	}
	REQUIRE(exceptionHappened);
}

TEST_CASE("CopyMatrixToBlock") {
	Matrix<float, 6, 7> mat(1);
	mat.block(1, 0, 4, 4) = Matrix<float, 4, 4>::zero();

	for (const auto& c : mat.block(1, 0, 4, 4))
		REQUIRE(0.f == c);
	bool exceptionHappened{ false };
	try {
		mat.block(1, 0, 4, 4) = Matrix<float, 3, 3>::zero();
	}
	catch (std::exception&) {
		exceptionHappened = true;
	}
	REQUIRE(exceptionHappened);
}

TEST_CASE("HadamardProduct") {
	Matrix<int, 2, 3> mat1{ 1, 2, 3, 4, 5, 6 };
	Matrix<int, 2, 3> mat2{ 23, -3, 4, 55, 622, 73 };
	auto prod = hadamard(mat1, mat2);
	for (auto it1 = mat1.begin(), it2 = mat2.begin(), it3 = prod.begin(); it1 != mat1.end(); ++it1, ++it2, ++it3) {
		REQUIRE((*it1) * (*it2) == *it3);
	}
}