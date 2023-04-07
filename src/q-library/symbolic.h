#pragma once
#include <string>
#include <iosfwd>
#include <complex>
#include <vector>
#include <iostream>
#include "matrix.h"

namespace Q {

	class Variable {
	public:
		//Variable(char c) { name = c; }
		Variable(const std::string& s) : name(s) {}
		friend std::ostream& operator<<(std::ostream& out, const Variable& v) { return out << v.name; }
		friend bool operator==(const Variable& a, const Variable& b) = default;
		//private:
		std::string name;
	};

	class Number {
	public:
		Number() = default;
		template<std::integral T> Number(T value) : value(value) {};
		template<std::floating_point T>	Number(T value) : value(value) {};
		template<std::floating_point T>	Number(const std::complex<T> value) : value(value) {};

		Number& operator+=(const Number& n) { value += n.value; return *this; }
		Number operator+(const Number& n) const { return Number{ *this } += n; }
		Number& operator*=(const Number& n) { value *= n.value; return *this; }
		Number operator*(const Number& n) const { return Number{ *this } *= n; }
		friend bool operator==(const Number& a, const Number& b) = default;

		friend std::ostream& operator<<(std::ostream& out, const Number& v) {
			if (v.value.imag() == 0.) return out << v.value.real();
			return out << '(' << v.value.real() << '+' << v.value.imag() << 'i' << ')';
		}
		double real() const { return value.real(); }
	private:
		std::complex<double> value;
	};



	class Sum;

	class Product {
	public:
		Product() = default;
		Product(const Number& n) { numbers.push_back(n); }

		Product& operator*=(const Number& v) { numbers.push_back(v); return *this; }
		Product& operator*=(const Variable& v) { variables.push_back(v); return *this; }
		Product& operator*=(const Sum& v);
		Product& operator*=(const Product& v) {
			numbers.insert(numbers.end(), v.numbers.begin(), v.numbers.end());
			variables.insert(variables.end(), v.variables.begin(), v.variables.end());
			sums.insert(sums.end(), v.sums.begin(), v.sums.end());
			return *this;
		}

		Product& simplify();

		// only yields correct result after call to simplify
		//bool isOneSimplified() const { return (numbers.size() == 0 || (numbers.size() == 1 && numbers[0] == 1.)) && variables.empty() && sums.empty(); }
		//bool isZeroSimplified() const { return numbers.size() == 1 && numbers[0] == 0 && variables.empty() && sums.empty(); }
		bool isNumeric() const { return variables.empty() && sums.empty(); }
		bool hasOnlySums() const { return numbers.empty() && variables.empty(); }
		Number numericProduct() const { return std::accumulate(numbers.begin(), numbers.end(), Number{ 1 }, std::multiplies{}); }
		size_t numTerms() const { return numbers.size() + variables.size() + sums.size(); }

		friend std::ostream& operator<<(std::ostream& out, const Product& v);

		//private:
		std::vector<Number> numbers;
		std::vector<Variable> variables;
		std::vector<Sum> sums;
	};

	class Sum {
	public:
		Sum() = default;
		Sum(const Number& n) { numbers.push_back(n); }

		Sum& operator+=(const Number& v) { numbers.push_back(v); return *this; }
		Sum& operator+=(const Variable& v) { variables.push_back(v); return *this; }
		Sum& operator+=(const Product& v);
		Sum& operator+=(const Sum& v) {
			numbers.insert(numbers.end(), v.numbers.begin(), v.numbers.end());
			variables.insert(variables.end(), v.variables.begin(), v.variables.end());
			products.insert(products.end(), v.products.begin(), v.products.end());
			return *this;
		}

		Sum& simplify();

		// only yields correct result after call to simplify
		//bool isZeroSimplified() const { return numbers.empty() && variables.empty() && products.empty(); }
		//bool isOneSimplified() const { return numbers.size() == 1 && numbers[0] == 1 && variables.empty() && products.empty(); }
		bool isNumeric() const { return variables.empty() && products.empty(); }
		bool hasOnlyProducts() const { return numbers.empty() && variables.empty(); }
		Number numericSum() const { return std::accumulate(numbers.begin(), numbers.end(), Number{}); }
		friend std::ostream& operator<<(std::ostream& out, const Sum& v);
		size_t numTerms() const { return numbers.size() + variables.size() + products.size(); }

		//private:
		std::vector<Number> numbers;
		std::vector<Variable> variables;
		std::vector<Product> products;
	};

	std::ostream& operator<<(std::ostream& out, const Product& v) {
		bool printedOne = false;
		if (!v.numbers.empty()) {
			out << v.numbers[0];
			for (size_t i = 1; i < v.numbers.size(); ++i) { out << '*' << v.numbers[i]; }
			printedOne = true;
		}
		if (!v.variables.empty()) {
			if (printedOne) out << '*';
			out << v.variables[0];
			for (size_t i = 1; i < v.variables.size(); ++i) { out << '*' << v.variables[i]; }
			printedOne = true;
		}
		if (!v.sums.empty()) {
			if (printedOne) out << '*';
			const bool needParentheses = (printedOne || v.sums.size() > 1) && v.sums[0].numTerms() > 0;
			if (needParentheses) out << '(';
			out << v.sums[0];
			if (needParentheses) out << ')';
			for (size_t i = 1; i < v.sums.size(); ++i) {
				if (v.sums[i].numTerms() > 1)
					out << "*(" << v.sums[i] << ')';
				else
					out << "*" << v.sums[i];
			}
			printedOne = true;
		}
		return out;
	}

	std::ostream& operator<<(std::ostream& out, const Sum& v) {
		bool printedOne = false;
		if (!v.numbers.empty()) {
			out << v.numbers[0];
			for (size_t i = 1; i < v.numbers.size(); ++i) { out << '+' << v.numbers[i]; }
			printedOne = true;
		}
		if (!v.variables.empty()) {
			if (printedOne) out << '+';
			out << v.variables[0];
			for (size_t i = 1; i < v.variables.size(); ++i) { out << '+' << v.variables[i]; }
			printedOne = true;
		}
		if (!v.products.empty()) {
			if (printedOne) out << '+';
			out << v.products[0];
			for (size_t i = 1; i < v.products.size(); ++i) { out << '+' << v.products[i]; }
			printedOne = true;
		}
		if (!printedOne) { out << "0"; }
		return out;
	}

	Product& Product::operator*=(const Sum& v) {
		if (v.numTerms() == 1) {
			for (const auto& c : v.numbers) (*this) *= c;
			for (const auto& c : v.variables) (*this) *= c;
			for (const auto& c : v.products) (*this) *= c;
		}
		else {
			sums.push_back(v);
		}
		return *this;
	}

	Sum& Sum::operator+=(const Product& v) {
		if (v.numTerms() == 1) {
			for (const auto& c : v.numbers) (*this) += c;
			for (const auto& c : v.variables) (*this) += c;
			for (const auto& c : v.sums) (*this) += c;
		}
		else {
			products.push_back(v);
		}
		return *this;
	}
	template <typename T>
	struct is_complex : std::false_type {};

	template <std::floating_point T>
	struct is_complex<std::complex<T>> : std::true_type {};

	class Term;

	template <typename T, typename... U>
	concept IsAnyOf = (std::same_as<T, U> || ...);
	template <typename T>
	concept Expression = IsAnyOf<T, Number, Variable, Product, Sum, Term>;
	template <typename T>
	concept Numeric = std::integral<T> || std::floating_point<T> || is_complex<T>::value;

	template<Expression T, Expression U>
	Sum operator+(const T& t, const U& u) {
		Sum sum;
		return (sum += t) += u;
	}

	template<Expression T, Numeric U>
	Sum operator+(const T& t, const U& u) { return t + Number{ u }; }

	template<Expression T, Numeric U>
	Sum operator+(const U& u, const T& t) { return t + Number{ u }; }

	template<Expression T, Expression U>
	Product operator*(const T& t, const U& u) {
		Product product;
		return (product *= t) *= u;
	}

	template<Expression T, Numeric U>
	Product operator*(const T& t, const U& u) { return t * Number{ u }; }

	template<Expression T, Numeric U>
	Product operator*(const U& u, const T& t) { return t * Number{ u }; }

	// Simplifiying cookbook
	//  1. All numbers are accumulated
	//  2. All sub-terms are simplified. All purely-numeric sub-terms are absorbed into this expression. 

	Product& Product::simplify() {
		Number n = 1;
		for (auto& sum : sums) {
			sum.simplify();
			if (sum.isNumeric()) {
				n *= sum.numericSum();
			}
			else if (sum.numTerms() == 1) {
				(*this) *= sum;
			}
		}
		n *= numericProduct();
		numbers.clear();
		if (n != 1)
			numbers.push_back(n);
		if (n == 0) {
			variables.clear();
			sums.clear();
		}
		std::erase_if(sums, [](const auto& sum) { return sum.isNumeric() || sum.numTerms() == 1; });
		return *this;
	}

	Sum& Sum::simplify() {
		Number n = 0;
		for (auto& product : products) {
			product.simplify();
			if (product.isNumeric()) {
				n += product.numericProduct();
			}
			else if (product.numTerms() == 1) {
				(*this) += product;
			}
		}
		n += numericSum();
		numbers.clear();
		if (n != 0)
			numbers.push_back(n);
		std::erase_if(products, [](const auto& product) { return product.isNumeric() || product.numTerms() == 1;  });
		return *this;
	}

	Sum& simplify(Sum& sum) { return sum.simplify(); }
	Product& simplify(Product& product) { return product.simplify(); }

	template<class Expr>
	Expr simplified(const Expr& expr) {
		auto e = expr;
		return simplify(e);
	}

	template<int m, int n>
	Math::Matrix<Term, m, n>& simplify(Math::Matrix<Term, m, n>& matrix) {
		for (auto& i : matrix) i.simplify();
		return matrix;
	}




	template<Numeric scalar, size_t m, size_t n, size_t p>
	constexpr Math::Matrix<Term, m, p> operator*(const Math::Matrix<scalar, m, n>& a, const Math::Matrix<Term, n, p>& b) {
		assert(a.cols() == b.rows());
		Math::Matrix<Term, m, p> result(a.shape() * b.shape());

		for (size_t i = 0; i < a.rows(); ++i) {
			for (size_t j = 0; j < b.cols(); ++j) {
				Term value{};
				for (size_t k = 0; k < a.cols(); ++k)
					value += a(i, k) * b(k, j);
				result(i, j) = value;
			}
		}
		return result;
	}


	template<Numeric scalar, size_t m, size_t n, size_t p>
	constexpr Math::Matrix<Term, m, p> operator*(const Math::Matrix<Term, m, n>& a, const Math::Matrix<scalar, n, p>& b) {
		assert(a.cols() == b.rows());
		Math::Matrix<Term, m, p> result(a.shape() * b.shape());

		for (size_t i = 0; i < a.rows(); ++i) {
			for (size_t j = 0; j < b.cols(); ++j) {
				Term value{};
				for (size_t k = 0; k < a.cols(); ++k)
					value += a(i, k) * b(k, j);
				result(i, j) = value;
			}
		}
		return result;
	}

	class Term : public Sum {
	public:
		Term() = default;
		template<Expression T>
		Term(const T& ex) { this->operator+=(ex); }
		template<Numeric T>
		Term(const T& ex) { this->operator+=(ex); }
	};


	// Sums have an advantage for Term-Term operations. In the case of mutiplying products, this can lead to
	// unnecessarily nested expressions. If both terms have exactly one summand, the product can be directly 
	// expressed as the product of the two summands. 
	Term operator*(const Term& t, const Term& u) {
		if (t.numTerms() == 1 && u.numTerms() == 1) {
			Product product;
			for (const auto& c : t.numbers) product *= c;
			for (const auto& c : t.variables) product *= c;
			for (const auto& c : t.products) product *= c;
			for (const auto& c : u.numbers) product *= c;
			for (const auto& c : u.variables) product *= c;
			for (const auto& c : u.products) product *= c;
			Term result;
			result += product;
			return result;
		}
		return Term{ static_cast<Sum>(t) * u }; // use normal operator overload for sums instead
	}




	template<int m, int n>
	Math::Matrix<Term, m, n> generateSymbolMatrix(const std::string& name) {
		Math::Matrix<Term, m, n> matrix;
		for (size_t i = 0; i < m; ++i) {
			for (size_t j = 0; j < n; ++j) {
				matrix(i, j) = Variable(name + std::to_string(i) + std::to_string(j));
			}
		}
		return matrix;
	}




	Math::Matrix<Term> generateSymbolMatrix(int m, int n, const std::string& name) {
		Math::Matrix<Term> matrix(m, n);
		for (size_t i = 0; i < m; ++i) {
			for (size_t j = 0; j < n; ++j) {
				matrix(i, j) = Variable(name + std::to_string(i) + std::to_string(j));
			}
		}
		return matrix;
	}

	template<int n>
	Math::Vector<Term, n> generateSymbolVector(const std::string& name) {
		Math::Vector<Term, n> vector;
		for (size_t j = 0; j < n; ++j) {
			vector[j] = Variable(name + std::to_string(j));
		}
		return vector;
	}

	constexpr Math::Matrix<Term> generateSymbolVector(int n, const std::string& name) {
		Math::Matrix<Term> vector(n, 1);
		for (size_t j = 0; j < n; ++j) {
			vector(j, 0) = Variable(name + std::to_string(j));
		}
		return vector;
	}

}
