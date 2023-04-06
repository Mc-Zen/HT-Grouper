#pragma once

namespace Q {


	/// @brief Binary number with modulo 2 arithmetic, 1+1=0.
	class Binary {
	public:
		constexpr Binary() = default;
		constexpr Binary(bool value) : value(value) {};
		explicit constexpr Binary(int value) : value(!!value) {};

		constexpr Binary& operator+=(const Binary& a) { value ^= a.value; return *this; }
		constexpr Binary& operator-=(const Binary& a) { return this->operator+=(a); }
		constexpr Binary& operator*=(const Binary& a) { value &= a.value; return *this; }
		constexpr Binary& operator|=(const Binary& a) { value |= a.value; return *this; }
		constexpr friend Binary operator+(const Binary& a, const Binary& b) { return Binary{ a } += b; }
		constexpr friend Binary operator-(const Binary& a, const Binary& b) { return Binary{ a } -= b; }
		constexpr friend Binary operator*(const Binary& a, const Binary& b) { return Binary{ a } *= b; }
		constexpr friend Binary operator|(const Binary& a, const Binary& b) { return Binary{ a } |= b; }

		constexpr operator int() const { return value; }
		constexpr int toInt() const { return value; }
		constexpr Binary& negate() { value = !value; return *this; }
		constexpr Binary operator~() const { return Binary{ !value }; }

	private:
		int value{};
	};
}
