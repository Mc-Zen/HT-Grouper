
#pragma once
#include <iterator>


namespace Q {

	int linearIndexFromNChoose2(int n, int i, int j) {
		return n * (n - 1) / 2 - (n - i) * (n - i - j) / 2 + j - i - 1;
	}

	std::pair<int, int> nChoose2FromLinearIndex(int n, int index) {
		const auto i = static_cast<int>(n - 2 - std::floor(std::sqrt(4 * n * (n - 1) - 7 - 8 * index) / 2.0 - 0.5));
		return { i, index + i + 1 - n * (n - 1) / 2 + (n - i) * (n - i - 1) / 2 };
	}

	class NChoose2Iterator {
	public:
		using iterator_concept = std::random_access_iterator_tag;
		using iterator_category = std::random_access_iterator_tag;
		using difference_type = std::ptrdiff_t;
		using size_type = size_t;
		using value_type = std::pair<int, int>;
		using pointer = value_type*;
		using reference = value_type&;

		constexpr NChoose2Iterator() noexcept = default;
		explicit NChoose2Iterator(int n, size_type index) noexcept : n(n), index(index) {
		}

		value_type operator*() const noexcept { return nChoose2FromLinearIndex(n, index); }

		constexpr NChoose2Iterator& operator++() noexcept { ++index; return *this; }
		constexpr NChoose2Iterator operator++(int) noexcept { NChoose2Iterator tmp = *this; ++index; return tmp; }
		constexpr NChoose2Iterator& operator--() noexcept { --index; return *this; }
		constexpr NChoose2Iterator operator--(int) noexcept { NChoose2Iterator tmp = *this; --index; return tmp; }
		constexpr NChoose2Iterator& operator+=(const difference_type offset) noexcept { index += offset; return *this; }
		constexpr NChoose2Iterator operator+(const difference_type offset) const noexcept { NChoose2Iterator tmp = *this; return tmp += offset; }
		constexpr NChoose2Iterator& operator-=(const difference_type offset) noexcept { index -= offset; return *this; }
		constexpr NChoose2Iterator operator-(const difference_type offset) const noexcept { NChoose2Iterator tmp = *this; return tmp -= offset; }
		constexpr difference_type operator-(const NChoose2Iterator& other) const noexcept { return index - other.index; }
		value_type operator[](const difference_type offset) const noexcept { return *(*this + offset); }
		constexpr std::strong_ordering operator<=>(const NChoose2Iterator& right) const noexcept { return index <=> right.index; }
		constexpr friend bool operator==(const NChoose2Iterator& a, const NChoose2Iterator& b) noexcept { return a.index == b.index && a.n == b.n; }
		constexpr friend NChoose2Iterator operator+(const difference_type offset, NChoose2Iterator a) noexcept { return a += offset; }

	private:
		int n{ 2 };
		int index{};
	};


	auto iterateThroughNChoose2(int n) {
		struct NChoose2 {
			constexpr NChoose2(int n) : n(n), last(n* (n - 1) / 2) {}
			auto begin() const { return  NChoose2Iterator{ n, 0 }; }
			auto end() const { return NChoose2Iterator{ n, last }; }
			auto cbegin() const { return begin(); }
			auto cend() const { return end(); }

			int n{};
			size_t last{};
		};

		return NChoose2{ n };
	}


}
