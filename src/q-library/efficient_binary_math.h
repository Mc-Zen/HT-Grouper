#pragma once
#include <cstdint>
#include "matrix.h"
#include "binary.h"

namespace Q::efficient {

	template<int n>
	class BinaryVector {
		static_assert(n < 65, "Only binary vectors with up to 64 entries are supported");
	public:
		constexpr BinaryVector() = default;
		constexpr BinaryVector(uint64_t value) : rep(value & mask) {}

		constexpr uint64_t operator()() const { return rep; }

		constexpr BinaryVector& operator+=(const BinaryVector& a) { rep ^= a.rep; return *this; }
		constexpr BinaryVector& operator*=(const BinaryVector& a) { rep &= a.rep; return *this; }
		constexpr BinaryVector& operator|=(const BinaryVector& a) { rep |= a.rep; return *this; }
		constexpr friend BinaryVector operator+(const BinaryVector& a, const BinaryVector& b) { return BinaryVector{ a } += b; }
		constexpr friend BinaryVector operator*(const BinaryVector& a, const BinaryVector& b) { return BinaryVector{ a } *= b; }
		constexpr friend BinaryVector operator|(const BinaryVector& a, const BinaryVector& b) { return BinaryVector{ a } |= b; }

		constexpr BinaryVector operator~() const { return { ~rep & mask }; }

		constexpr uint64_t dot(const BinaryVector& a) const {
			auto u = std::popcount(rep & a.rep) & 1;
			//println("{},{},{}, {}", rep, a.rep, rep & a.rep, u);
			return u;
		}
		constexpr uint64_t bitCount() const { return std::popcount(rep); }
		friend constexpr bool operator==(BinaryVector a, BinaryVector b) = default;
		constexpr uint64_t get(size_t i) const { assert(i >= 0 && i < n && "Invaild index"); return (rep & (1ULL << i)) != 0; }
		constexpr void set(size_t i, uint64_t value) { assert(i >= 0 && i < n && "Invaild index"); rep |= (static_cast<uint64_t>(value) << i); }



		auto toVector() const {
			Math::Vector<Q::Binary, n> result;
			auto r = rep;
			for (size_t i = 0; i < n; ++i) {
				result[i] = r & 1;
				r >>= 1;
			}
			return result;
		}

	private:
		static constexpr uint64_t mask = (1ULL << n) - 1;
		uint64_t rep{};
	};


	template<int n>
	auto createBinaryVector(BinaryVector<n> bitstring) {
		auto rep = bitstring();
		Math::Vector<Q::Binary, n> result;
		for (size_t i = 0; i < n; ++i) {
			result[i] = rep & 1;
			rep >>= 1;
		}
		return result;
	}

	template<int n>
	auto toBitstringInteger(const Math::Vector<Q::Binary, n>& vec) {
		BinaryVector<n> output{};
		for (size_t i = 0; i < n; ++i) {
			output |= static_cast<uint64_t>(vec[i]) << i;
		}
		return output;
	}


	template<int m, int n>
	struct BinaryColMatrix {
		using Column = BinaryVector<m>;
		std::array<Column, n> cols{};

		static_assert(m < 4 * sizeof(Column));

		constexpr Column operator[](int col) const { return cols[col]; }
		constexpr Column& operator[](int col) { return cols[col]; }

		constexpr BinaryColMatrix() = default;
		explicit constexpr BinaryColMatrix(const Math::Matrix<Q::Binary, m, n>& mat) {
			for (size_t col = 0; col < n; ++col) {
				Column colVec{};
				for (size_t row = 0; row < m; ++row) {
					colVec.set(row, static_cast<uint64_t>(mat(row, col)));
				}
				cols[col] = colVec;
			}
		}

		constexpr auto toMatrix() const {
			Math::Matrix<Q::Binary, m, n> mat;
			for (size_t row = 0; row < m; ++row) {
				for (size_t col = 0; col < n; ++col) {
					mat(row, col) = cols[col]() & (1ULL << row);
				}
			}
			return mat;
		}
	};

	template<int m, int n>
	constexpr efficient::BinaryVector<n> operator*(efficient::BinaryVector<m> vector, const BinaryColMatrix<m, n>& matrix) {
		efficient::BinaryVector<n> result{};
		for (size_t col = 0; col < n; ++col) {
			result.set(col, matrix[col].dot(vector));
		}
		return result;
	}

	template<int m, int n>
	struct BinaryRowMatrix {
		using Row = BinaryVector<n>;
		std::array<Row, m> rows{};

		static_assert(m < 4 * sizeof(Row));


		constexpr BinaryRowMatrix() = default;

		explicit constexpr BinaryRowMatrix(const Math::Matrix<Binary, m, n>& mat) {
			for (size_t row = 0; row < m; ++row) {
				Row rowVec{};
				for (size_t col = 0; col < n; ++col) {
					rowVec.set(col, static_cast<uint64_t>(mat(row, col)));
				}
				rows[row] = rowVec;
			}
		}

		constexpr Row operator[](int row) const { return rows[row]; }
		constexpr Row& operator[](int row) { return rows[row]; }


		constexpr auto toMatrix() const {
			Math::Matrix<Binary, m, n> mat;
			for (size_t row = 0; row < m; ++row) {
				for (size_t col = 0; col < n; ++col) {
					mat(row, col) = rows[row].get(col);
				}
			}
			return mat;
		}

	};


	template<int m, int n>
	constexpr efficient::BinaryVector<m> operator*(const BinaryRowMatrix<m, n>& matrix, efficient::BinaryVector<n> vector) {
		BinaryVector<m> result{};
		for (size_t row = 0; row < m; ++row) {
			result.set(row, matrix[row].dot(vector));
		}
		return result;
	}
}
