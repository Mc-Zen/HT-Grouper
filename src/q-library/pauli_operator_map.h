
#pragma once
#include "binary_pauli.h"
#include <vector>

namespace Q {



	/// @brief Compact binary representation of an n-qubit pauli operator. For each qubit
	///        two bits are used to represent a pauli operator of {I,X,Z,Y}
	///           00 -> I,    01 -> X,    10 -> Z,    11 -> Y
	///        Nice feature: taking them as indices, this allows a dense representation of the 
	///        n-qubit Pauli group. 
	struct PauliIndex {
		constexpr PauliIndex() = default;

		explicit constexpr PauliIndex(int numQubits, uint64_t index) : numQubits(numQubits), index(index) {}

		explicit constexpr PauliIndex(const std::string_view& str) {
			assert(str.size() == numQubits);
			for (int i = 0; i < numQubits; ++i) {
				index <<= 2;
				switch (str[i]) {
				case 'X': index |= 0b01; break;
				case 'Y': index |= 0b11; break;
				case 'Z': index |= 0b10; break;
				default: break;
				}
			}
		}

		explicit constexpr PauliIndex(const Pauli& pauli) {
			assert(pauli.numQubits() == numQubits);
			for (int i = 0; i < numQubits; ++i) {
				index <<= 1;
				index |= pauli.z(i);
				index <<= 1;
				index |= pauli.x(i);
			}
		}

		std::string toString() const {
			static constexpr std::array<char, 4> c{ 'I','X','Z','Y' };
			std::string result;
			for (int i = 0; i < numQubits; ++i) {
				result += c[(index >> (2 * (numQubits - i - 1))) & 0b11];
			}
			return result;
		}


		int numQubits{};
		uint64_t index{};
	};


	template<class T>
	class PauliOperatorMap {
	public:


		template<class U>
		class PauliIndexIterator {
		public:
			using iterator_concept = std::contiguous_iterator_tag;
			using iterator_category = std::random_access_iterator_tag;
			using difference_type = std::ptrdiff_t;
			using size_type = size_t;
			using value_type = std::pair<PauliIndex, std::remove_const_t<U>>;
			using pointer = U*;
			using reference = std::pair<PauliIndex, U&>;

			constexpr PauliIndexIterator() noexcept = default;
			constexpr explicit PauliIndexIterator(int numQubits, pointer ptr, size_type offset = 0) noexcept
				: numQubits(numQubits), ptr(ptr + offset), index(offset) {}

			constexpr reference operator*() const noexcept { return { PauliIndex{ numQubits, index }, *ptr }; }
			//constexpr pointer operator->() const noexcept { return ptr; }
			constexpr PauliIndexIterator& operator++() noexcept { ++ptr; ++index; return *this; }
			constexpr PauliIndexIterator operator++(int) noexcept { PauliIndexIterator tmp = *this; ++ptr; ++index; return tmp; }
			constexpr PauliIndexIterator& operator--() noexcept { --ptr; --index; return *this; }
			constexpr PauliIndexIterator operator--(int) noexcept { PauliIndexIterator tmp = *this; --ptr; --index; return tmp; }
			constexpr PauliIndexIterator& operator+=(const difference_type offset) noexcept { ptr += offset; index += offset;  return *this; }
			constexpr PauliIndexIterator operator+(const difference_type offset) const noexcept { PauliIndexIterator tmp = *this; return tmp += offset; }
			constexpr PauliIndexIterator& operator-=(const difference_type offset) noexcept { ptr -= offset; index -= offset; return *this; }
			constexpr PauliIndexIterator operator-(const difference_type offset) const noexcept { PauliIndexIterator tmp = *this; return tmp -= offset; }
			constexpr difference_type operator-(const PauliIndexIterator& other) const noexcept { return ptr - other.ptr; }
			constexpr reference operator[](const difference_type offset) const noexcept { return ptr[offset]; }
			constexpr std::strong_ordering operator<=>(const PauliIndexIterator& right) const noexcept { return ptr <=> right.ptr; }
			constexpr friend bool operator==(const PauliIndexIterator& a, const PauliIndexIterator& b) noexcept { return a.ptr == b.ptr; }
			constexpr friend PauliIndexIterator operator+(const difference_type offset, const PauliIndexIterator a) noexcept { return a += offset; }

		private:
			int numQubits{};
			pointer ptr{};
			size_type index{};
		};

		template<class U>
		struct PauliEnumerator {
			U& map;

			auto begin() const { return PauliIndexIterator<const T>{map.map.data(), 0}; }
			auto end() const { return PauliIndexIterator<const T>{map.map.data(), map.map.size()}; }
			auto cbegin() const { return begin(); }
			auto cend() const { return end(); }
			auto begin() { return PauliIndexIterator<T>{map.map.data(), 0}; }
			auto end() { return PauliIndexIterator<T>{map.map.data(), map.map.size()}; }
		};

		using value_type = T;
		using size_type = size_t;
		using difference_type = std::ptrdiff_t;

		using iterator = typename std::vector<value_type>::iterator;
		using const_iterator = typename std::vector<value_type>::const_iterator;
		using reverse_iterator = typename std::vector<value_type>::reverse_iterator;
		using const_reverse_iterator = typename std::vector<value_type>::const_reverse_iterator;


		explicit constexpr PauliOperatorMap(int numQubits) : numQubits(numQubits), map(pow4(numQubits)) {}
		explicit constexpr PauliOperatorMap(int numQubits, const T& value) : numQubits(numQubits), map(pow4(numQubits), value) {}

		constexpr T& operator[](const std::string_view& s) { return map[toIndex(s)]; }
		constexpr const T& operator[](const std::string_view& s) const { return map[toIndex(s)]; }
		constexpr T& operator[](const PauliIndex& s) {
			assert(s == numQubits);
			return map[s.index];
		}
		constexpr const T& operator[](const PauliIndex& s) const {
			assert(s == numQubits);
			return map[s.index];
		}

		constexpr iterator begin() { return map.begin(); }
		constexpr iterator end() { return map.end(); }
		constexpr const_iterator begin() const { return map.begin(); }
		constexpr const_iterator end() const { return map.end(); }
		constexpr const_iterator cbegin() const { return begin(); }
		constexpr const_iterator cend() const { return end(); }

		constexpr reverse_iterator rbegin() { return map.rbegin(); }
		constexpr reverse_iterator rend() { return map.rend(); }
		constexpr const_reverse_iterator rbegin() const { return map.rbegin(); }
		constexpr const_reverse_iterator rend() const { return map.rend(); }
		constexpr const_reverse_iterator crbegin() const { return rbegin(); }
		constexpr const_reverse_iterator crend() const { return rend(); }

		// Allows to enumerate entries using
		//    for (auto&& [pauliIndex, value] : map.enumerate()) {
		//    }
		// which returns pairs of a PauliIndex and the value which is contained for this index. 
		auto enumerate() const { return PauliEnumerator{ *this }; }
		auto enumerate() { return PauliEnumerator{ *this }; }


	private:
		size_t toIndex(const std::string_view& s) const {
			assert(s.size() == numQubits);
			return PauliIndex{ numQubits, s }.index;
		}

		std::vector<T> map;
		int numQubits{};
	};

}
