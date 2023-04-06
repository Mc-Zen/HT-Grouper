#pragma once
#include <complex>
#include <numbers>
#include "matrix.h"
#include "special_math.h"


namespace Q {

	using real = double;
	using scalar = std::complex<real>;

	using namespace std::literals::complex_literals;
	inline constexpr scalar sqrt2 = 1.4142135623730950488016887242097;
	inline constexpr scalar invsqrt2 = 1. / sqrt2;


	template<int numQubits>
	using Operator = Math::Matrix<scalar, pow2(numQubits), pow2(numQubits)>;

	/// @brief 1-Qubit Operator
	using Op1 = Operator<1>;
	/// @brief 2-Qubit Operator
	using Op2 = Operator<2>;
	/// @brief 3-Qubit Operator
	using Op3 = Operator<3>;


	template<class Operator>
	struct OperatorTraits { static constexpr int numQubits() { return log2OfPowerOf2(static_cast<int>(Operator::rows())); } };

	//template<> struct OperatorTraits<Operator<1>> { static constexpr int numQubits() { return 1; } };
	//template<> struct OperatorTraits<Operator<1>> { static constexpr int numQubits() { return 1; } };
	//template<> struct OperatorTraits<Operator<2>> { static constexpr int numQubits() { return 2; } };
	//template<> struct OperatorTraits<Operator<3>> { static constexpr int numQubits() { return 3; } };
	//template<> struct OperatorTraits<Operator<4>> { static constexpr int numQubits() { return 4; } };




	namespace Gates {
		inline constexpr Op1 I{ 1,0,0,1 };
		inline constexpr Op1 H{ invsqrt2,invsqrt2,invsqrt2,-invsqrt2 };
		inline constexpr Op1 X{ 0,1,1,0 };
		inline constexpr Op1 Y{ 0,-1i,1i,0 };
		inline constexpr Op1 Z{ 1,0,0,-1 };
		inline constexpr Op1 S{ 1,0,0,1i };
		inline constexpr Op1 SDG{ 1,0,0,-1i };
		inline constexpr Op1 T{ 1,0,0,sqrt2 + 1i * sqrt2 };
		inline constexpr Op2 CX{ 1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0 };
		inline constexpr Op2 CZ{ 1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,-1 };

		inline Op1 U(real theta, real phi, real lambda) {
			const auto sinTheta2 = std::sin(0.5 * theta);
			const auto cosTheta2 = std::cos(0.5 * theta);
			return Op1{
				cosTheta2, -std::polar<real>(1, lambda) * sinTheta2,
				std::polar<real>(1,phi) * sinTheta2, std::polar<real>(1, phi + lambda) * cosTheta2
			};
		}
	}

	//constexpr Op2 tensor(const Op1& op1, const Op1& op2) {
	//	Op2 op{};
	//	op.block(0, 0, 2, 2) = op2 * op1(0, 0);
	//	op.block(0, 2, 2, 2) = op2 * op1(0, 1);
	//	op.block(2, 0, 2, 2) = op2 * op1(1, 0);
	//	op.block(2, 2, 2, 2) = op2 * op1(1, 1);
	//	return op;
	//}

	/// @brief Tensor product of two operators A \otimes B
	template<class Op1, class Op2>
	constexpr auto tensorProduct(const Op1& A, const Op2& B) {
		Operator<OperatorTraits<Op1>::numQubits() + OperatorTraits<Op2>::numQubits()> result{};
		constexpr auto n = B.cols();
		for (size_t i = 0; i < A.cols(); ++i) {
			for (size_t j = 0; j < A.rows(); ++j) {
				result.block(i * n, j * n, n, n) = A(i, j) * B;
			}
		}
		return result;
	}

	template<class Op1, class Op2>
	constexpr auto operator%(const Op1& A, const Op2& B) {
		Operator<OperatorTraits<Op1>::numQubits() + OperatorTraits<Op2>::numQubits()> result{};
		constexpr auto n = B.cols();
		for (size_t i = 0; i < A.cols(); ++i) {
			for (size_t j = 0; j < A.rows(); ++j) {
				result.block(i * n, j * n, n, n) = A(i, j) * B;
			}
		}
		return result;
	}

	template<int numQubits>
	auto tensorProduct(const std::array<Op1, numQubits>& operators) {
		static constexpr auto sum = [](auto&&... items) { return (... % items); };
		static constexpr auto accumulate2 = [](auto&& fn, const auto& u) {
			return std::apply([fn = std::forward<decltype(fn)>(fn)](auto... item) { return fn(item...); }, u);
		};
		return accumulate2(sum, operators);
	}


	/// @brief Generates a controlled gate from the given operator
	constexpr Op2 controlled(const Op1& op) {
		Op2 result{};
		result.block(0, 0, 2, 2) = Gates::I;
		result.block(2, 2, 2, 2) = op;
		return result;
	}



	template<class T, int m, int n>
	constexpr Math::Matrix<T, n, m> dagger(const Math::Matrix<T, m, n>& mat) {
		auto z = mat.transpose();
		for (auto& c : z) {
			c = std::conj(c);
		}
		return z;
	}


	/// @brief Compute how a given operator transforms under the transformation of another operator. 
	template<class Op>
	constexpr Op transformOperator(const Op& op, const Op& transform) {
		return transform * op * dagger(transform);
	}

	/// @brief Compute how a given operator transforms under the transformation of another operator. 
	template<class Op>
	constexpr Op commutator(const Op& op1, const Op& op2) {
		return op1 * op2 - op2 * op1;
	}


	/// @brief Very rough factorizsation of a tensor product of two pauli matrices. Also works with
	///        Similar operators like H, S, T
	/// @param op 2-Qubit gate
	/// @return two 1-Qubit gates, so that tensorProduct(result.first, result.second) == op
	std::pair<Op1, Op1> pauliFactorize(const Op2& op) {
		Op1 op1{ (op(0,0) == 0. && op(0,1) == 0.) ? op.block(0,2,2,2) : op.block(0,0,2,2) };
		Op1 op2;
		if (op1(0, 0) != 0.) {
			op2(0, 0) = op(0, 0) / op1(0, 0);
			op2(1, 0) = op(2, 0) / op1(0, 0);
			op2(0, 1) = op(0, 2) / op1(0, 0);
			op2(1, 1) = op(2, 2) / op1(0, 0);
		}
		else {
			op2(0, 0) = op(0, 1) / op1(0, 1);
			op2(1, 0) = op(2, 1) / op1(0, 1);
			op2(0, 1) = op(0, 3) / op1(0, 1);
			op2(1, 1) = op(2, 3) / op1(0, 1);
		}

		return { op2,op1 };
	}
}
