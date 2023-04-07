
#pragma once
#include <string>
#include <array>


namespace Q {

	/// @brief Complex binary phase in the form i^q using an exponent q = 0,1,2,3
	class BinaryPhase {
	public:
		constexpr BinaryPhase() = default;
		explicit(false) constexpr BinaryPhase(int n) { set(n); }

		constexpr BinaryPhase& operator+=(int n) { set(phase + n); return *this; }
		constexpr BinaryPhase& operator-=(int n) { return *this += (n * 3); }
		constexpr BinaryPhase& operator+=(BinaryPhase other) { return *this += other.phase; }
		constexpr BinaryPhase& operator-=(BinaryPhase other) { return *this -= other.phase; }
		constexpr friend BinaryPhase operator+(BinaryPhase a, BinaryPhase b) { return a += b; }
		constexpr friend BinaryPhase operator-(BinaryPhase a, BinaryPhase b) { return a -= b; }

		constexpr BinaryPhase& operator++() { return *this += 1; }
		constexpr BinaryPhase& operator--() { return *this += 3; }
		constexpr unsigned int operator()() const { return phase; }
		constexpr unsigned int toInt() const { return phase; }

		std::string toString() const { static constexpr std::array<std::string_view, 4> c{ "+","i","-","-i" };  return std::string(c[phase]); }

		constexpr bool isPlusMinus() const { return phase == 0 || phase == 2; }
		constexpr friend bool operator==(const BinaryPhase& a, const BinaryPhase& b) = default;

	private:
		constexpr void set(int n) { phase = n & 0b11; }
		unsigned int phase{};
	};

}


