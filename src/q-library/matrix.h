#pragma once

#include <array>
#include <cmath>
#include <numeric>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <iosfwd>
#include <cassert>
#include <vector>



#ifdef MATRIX_EXCEPTIONS
#define MATRIX_VERIFY(expr, msg, exception) if(!(expr)) throw exception(msg)
#else
#define MATRIX_VERIFY(expr, msg, exception) assert(msg && expr)
#endif

namespace Math {

	using Index = size_t;

	// concepts for identifying dimension combinations as vector dimensions (of special size)
	template<Index m, Index n>
	concept vector_dimensions = m == 1 || n == 1;

	template<Index m, Index n, Index dmin, Index dmax>
	concept vector_dimension_limits = m == 1 && dmin <= n && n <= dmax;

	template<Index m, Index n>
	concept vector_dimensions_1D = vector_dimension_limits<m, n, 1, 4> || vector_dimension_limits<n, m, 1, 4>;
	template<Index m, Index n>
	concept vector_dimensions_2D = vector_dimension_limits<m, n, 2, 4> || vector_dimension_limits<n, m, 2, 4>;
	template<Index m, Index n>
	concept vector_dimensions_3D = vector_dimension_limits<m, n, 3, 4> || vector_dimension_limits<n, m, 3, 4>;
	template<Index m, Index n>
	concept vector_dimensions_4D = vector_dimension_limits<m, n, 4, 4> || vector_dimension_limits<n, m, 4, 4>;


	struct Matrix_view_mismatch : public std::logic_error { using std::logic_error::logic_error; };
	struct Matrix_block_domain_error : public std::domain_error { using std::domain_error::domain_error; };
	struct Matrix_shape_error : public std::length_error { using std::length_error::length_error; };

	inline constexpr Index dynamic = std::numeric_limits<Index>::max();

	template<class T, Index m, Index n> class Matrix;

	template<class U>
	class Matrix_iterator {
	public:
		using iterator_concept = std::contiguous_iterator_tag;
		using iterator_category = std::random_access_iterator_tag;
		using difference_type = std::ptrdiff_t;
		using size_type = size_t;
		using value_type = std::remove_const_t<U>;
		using pointer = U*;
		using reference = U&;

		constexpr Matrix_iterator() noexcept = default;
		constexpr explicit Matrix_iterator(pointer ptr, size_type offset = 0) noexcept : ptr(ptr + offset) {}

		constexpr reference operator*() const noexcept { return *ptr; }
		constexpr pointer operator->() const noexcept { return ptr; }
		constexpr Matrix_iterator& operator++() noexcept { ++ptr; return *this; }
		constexpr Matrix_iterator operator++(int) noexcept { Matrix_iterator tmp = *this; ++ptr; return tmp; }
		constexpr Matrix_iterator& operator--() noexcept { --ptr; return *this; }
		constexpr Matrix_iterator operator--(int) noexcept { Matrix_iterator tmp = *this; --ptr; return tmp; }
		constexpr Matrix_iterator& operator+=(const difference_type offset) noexcept { ptr += offset; return *this; }
		constexpr Matrix_iterator operator+(const difference_type offset) const noexcept { Matrix_iterator tmp = *this; return tmp += offset; }
		constexpr Matrix_iterator& operator-=(const difference_type offset) noexcept { ptr -= offset; return *this; }
		constexpr Matrix_iterator operator-(const difference_type offset) const noexcept { Matrix_iterator tmp = *this; return tmp -= offset; }
		constexpr difference_type operator-(const Matrix_iterator& other) const noexcept { return ptr - other.ptr; }
		constexpr reference operator[](const difference_type offset) const noexcept { return ptr[offset]; }
		constexpr std::strong_ordering operator<=>(const Matrix_iterator& right) const noexcept { return ptr <=> right.ptr; }
		constexpr friend bool operator==(const Matrix_iterator& a, const Matrix_iterator& b) noexcept = default;
		constexpr friend Matrix_iterator operator+(const difference_type offset, Matrix_iterator a) noexcept { return a += offset; }

	private:
		pointer ptr{};
	};


	template<class U>
	class Matrix_stride_iterator {
	public:
		using iterator_concept = std::random_access_iterator_tag;
		using iterator_category = std::random_access_iterator_tag;
		using difference_type = std::ptrdiff_t;
		using size_type = size_t;
		using value_type = std::remove_const_t<U>;
		using pointer = U*;
		using reference = U&;

		constexpr Matrix_stride_iterator() noexcept = default;
		constexpr explicit Matrix_stride_iterator(pointer ptr, size_type stride = 1, size_type offset = 0) noexcept : ptr(ptr + offset), stride(stride) {}

		constexpr reference operator*() const noexcept { return *ptr; }
		constexpr pointer operator->() const noexcept { return ptr; }
		constexpr Matrix_stride_iterator& operator++() noexcept { ptr += stride; return *this; }
		constexpr Matrix_stride_iterator operator++(int) noexcept { Matrix_stride_iterator tmp = *this; ptr += stride; return tmp; }
		constexpr Matrix_stride_iterator& operator--() noexcept { ptr -= stride; return *this; }
		constexpr Matrix_stride_iterator operator--(int) noexcept { Matrix_stride_iterator tmp = *this; ptr -= stride; return tmp; }
		constexpr Matrix_stride_iterator& operator+=(const difference_type offset) noexcept { ptr += offset * stride; return *this; }
		constexpr Matrix_stride_iterator operator+(const difference_type offset) const noexcept { Matrix_stride_iterator tmp = *this; return tmp += offset * stride; }
		constexpr Matrix_stride_iterator& operator-=(const difference_type offset) noexcept { ptr -= offset * stride; return *this; }
		constexpr Matrix_stride_iterator operator-(const difference_type offset) const noexcept { Matrix_stride_iterator tmp = *this; return tmp -= offset * stride; }
		constexpr difference_type operator-(const Matrix_stride_iterator& other) const noexcept { return ptr - other.ptr; }
		constexpr reference operator[](const difference_type offset) const noexcept { return ptr[offset * stride]; }
		constexpr std::strong_ordering operator<=>(const Matrix_stride_iterator& right) const noexcept { return ptr <=> right.ptr; }
		constexpr friend bool operator==(const Matrix_stride_iterator& a, const Matrix_stride_iterator& b) noexcept { return a.ptr == b.ptr; }
		constexpr friend Matrix_stride_iterator operator+(const difference_type offset, const Matrix_stride_iterator a) noexcept { return a += offset; }

	private:
		pointer ptr{};
		size_type stride{ 1 };
	};


	template<class U>
	class Matrix_block_iterator {
	public:
		using iterator_concept = std::bidirectional_iterator_tag;
		using iterator_category = std::bidirectional_iterator_tag;
		using difference_type = std::ptrdiff_t;
		using size_type = size_t;
		using value_type = std::remove_const_t<U>;
		using pointer = U*;
		using reference = U&;

		constexpr Matrix_block_iterator() noexcept = default;
		constexpr explicit Matrix_block_iterator(pointer ptr, size_type row_length = 1, size_type jump = 1, size_type row_index = 0) noexcept
			: ptr(ptr), row_length(row_length), jump(jump), row_index(row_index) {}

		constexpr reference operator*() const noexcept { return *ptr; }
		constexpr pointer operator->() const noexcept { return ptr; }
		constexpr Matrix_block_iterator& operator++() noexcept {
			if (++row_index >= row_length) { ptr += jump; row_index = 0; }
			++ptr; return *this;
		}
		constexpr Matrix_block_iterator operator++(int) noexcept { Matrix_block_iterator tmp = *this; (*this)++; return tmp; }
		constexpr Matrix_block_iterator& operator--() noexcept {
			if (row_index-- == 0) { ptr -= jump; row_index = row_length - 1; }
			ptr--; return *this;
		}
		constexpr Matrix_block_iterator operator--(int) noexcept { Matrix_block_iterator tmp = *this; (*this)--; return tmp; }
		constexpr std::strong_ordering operator<=>(const Matrix_block_iterator& right) const noexcept { return ptr <=> right.ptr; }
		constexpr friend bool operator==(const Matrix_block_iterator& a, const Matrix_block_iterator& b) noexcept { return a.ptr == b.ptr; }

	private:
		pointer ptr{};
		size_type row_length{ 1 };
		size_type row_index{};
		size_type jump{ 1 };
	};


	template<class T>
	class Matrix_view {
	public:
		using value_type = T;
		using pointer = T*;
		using reference = T&;
		using size_type = size_t;

		using iterator = Matrix_block_iterator<T>;
		using const_iterator = Matrix_block_iterator<const T>;

		using reverse_iterator = std::reverse_iterator<iterator>;
		using const_reverse_iterator = std::reverse_iterator<const_iterator>;


		constexpr Matrix_view() noexcept = default;
		constexpr Matrix_view(T* ptr, size_type m, size_type n, size_type rows, size_type cols) noexcept
			: ptr(ptr), m(m), n(n), rows_(rows), cols_(cols) {
		}

		constexpr Matrix_view(const Matrix_view<T>&) noexcept = default;
		constexpr ~Matrix_view() noexcept = default;

		constexpr Matrix_view& operator=(const Matrix_view& mat) requires (!std::is_const_v<T>) {
			MATRIX_VERIFY(rows() == mat.rows() && cols() == mat.cols(), "Matrix view mismatch at Matrix_view::operator=(const Matrix_view&). Dimensions of target and destination need to match", Matrix_view_mismatch);
			std::copy(mat.begin(), mat.end(), begin());
			return *this;
		}

		//constexpr Matrix_view& operator=(const Matrix_view<const T>& mat) requires (!std::is_const_v<T>) {
		//	MATRIX_VERIFY(rows() == mat.rows() && cols() == mat.cols(), "Matrix view mismatch at Matrix_view::operator=(const Matrix_view&). Dimensions of target and destination need to match", Matrix_view_mismatch);
		//	std::copy(mat.begin(), mat.end(), begin());
		//	return *this;
		//}

		template<size_type m, size_type n>
		constexpr Matrix_view& operator=(const Matrix<T, m, n>& mat) requires (!std::is_const_v<T>) {
			MATRIX_VERIFY(rows() == mat.rows() && cols() == mat.cols(), "Matrix view mismatch at Matrix_view::operator=(const Matrix&). Dimensions of target and destination need to match", Matrix_view_mismatch);
			std::copy(mat.begin(), mat.end(), begin());
			return *this;
		}

		constexpr iterator begin() noexcept requires (!std::is_const_v<T>) { return iterator(ptr, cols(), n - cols()); }
		constexpr iterator end() noexcept requires (!std::is_const_v<T>) { return iterator(ptr + rows() * n, cols(), n - cols()); }
		constexpr const_iterator begin() const noexcept { return const_iterator(ptr, cols(), n - cols()); }
		constexpr const_iterator end() const noexcept { return const_iterator(ptr + rows() * n, cols(), n - cols()); }

		constexpr reverse_iterator rbegin() noexcept requires (!std::is_const_v<T>) { return reverse_iterator(end()); }
		constexpr reverse_iterator rend()  requires (!std::is_const_v<T>) { return reverse_iterator(begin()); }
		constexpr const_reverse_iterator rbegin() const noexcept { return const_reverse_iterator(end()); }
		constexpr const_reverse_iterator rend() const noexcept { return const_reverse_iterator(begin()); }


		constexpr size_type rows() const noexcept { return rows_; }
		constexpr size_type cols() const noexcept { return cols_; }

	private:
		T* ptr{}; // first element of this matrix view block
		size_type m{};
		size_type n{}; // dimensions of matrix that this view is pointing to 
		size_type rows_{};
		size_type cols_{}; // dimensions of this matrix view block

	};


	template<bool empty = false>
	struct Shape {
		Index m{};
		Index n{};
		constexpr void swap() { std::swap(m, n); }
		constexpr Shape swapped() const { auto copy = *this; copy.swap(); return copy; }
		friend constexpr Shape operator*(Shape d1, Shape d2) { return Shape{ d1.m, d2.n }; }
		friend constexpr bool operator==(const Shape& a, const Shape& b) = default;
	};

	template<>
	struct Shape<true> {
		constexpr void swap() const { /**/ }
		constexpr Shape swapped() const { return *this; }
		friend constexpr Shape operator*(Shape, Shape) { return {}; }
		friend constexpr bool operator==(const Shape& a, const Shape& b) = default;
	};


	template<class T, Index m = dynamic, Index n = dynamic>
	class Matrix {
	public:

		static constexpr Index m_ = m;
		static constexpr Index n_ = n;

		using value_type = T;
		using reference = T&;
		using const_reference = const T&;
		using size_type = Index;
		using difference_type = ptrdiff_t;

		//using storage_type = std::array<T, m* n>;
		using storage_type = std::conditional_t<m == dynamic || n == dynamic, std::vector<T>, std::array<T, m* n>>;
		static_assert(m > 0 || m == dynamic, "Row number m needs to be greater than zero");
		static_assert(n > 0 || m == dynamic, "Column number m needs to be greater than zero");
		static_assert((m == dynamic) == (n == dynamic), "Row and column number cannot be independantly dynamic");
		static constexpr bool is_dynamic = (m == dynamic);

		using iterator = Matrix_iterator<T>;
		using const_iterator = Matrix_iterator<const T>;
		using reverse_iterator = std::reverse_iterator<iterator>;
		using const_reverse_iterator = std::reverse_iterator<const_iterator>;

		using row_iterator = Matrix_iterator<T>;
		using const_row_iterator = Matrix_iterator<const T>;
		using col_iterator = Matrix_stride_iterator<T>;
		using const_col_iterator = Matrix_stride_iterator<const T>;


		//
		// Constructors
		//


		constexpr Matrix() = default;

		explicit constexpr Matrix(const T& value) requires (!is_dynamic) { fill(value); }

		explicit(false) constexpr Matrix(const std::initializer_list<T> elems) requires (!is_dynamic) {
			size_type h = std::min(elems.size(), size());
			std::copy(std::begin(elems), std::begin(elems) + h, begin());
			std::fill(begin() + h, end(), T{});
		}

		explicit constexpr Matrix(const storage_type& elems) requires (!is_dynamic) { std::copy(std::begin(elems), std::end(elems), begin()); }

		explicit constexpr Matrix(storage_type&& elems) noexcept requires (!is_dynamic) : data_{ std::move(elems) } {}

		explicit(false) constexpr Matrix(const Matrix_view<T>& mat_view) {
			if constexpr (is_dynamic) {
				resize(mat_view.rows(), mat_view.cols());
			}
			else {
				MATRIX_VERIFY(mat_view.rows() == rows() && mat_view.cols() == cols(), "Dimensions of block and matrix do not match in Matrix::operator=(const Matrix_view&)", Matrix_block_domain_error);
			}
			std::copy(mat_view.begin(), mat_view.end(), begin());
		}

		explicit(false) constexpr Matrix(const Matrix_view<const T>& mat_view) {
			if constexpr (is_dynamic) {
				resize(mat_view.rows(), mat_view.cols());
			}
			else {
				MATRIX_VERIFY(mat_view.rows() == rows() && mat_view.cols() == cols(), "Dimensions of block and matrix do not match in Matrix::operator=(const Matrix_view&)", Matrix_block_domain_error);
			}
			std::copy(mat_view.begin(), mat_view.end(), begin());
		}

		template<bool a>
		explicit constexpr Matrix(Shape<a> shape_) : shape_(shape_) {
			if constexpr (is_dynamic) {
				data_.resize(size());
			}
		}

		constexpr Matrix(Index m_, Index n_) requires is_dynamic : data_(m_* n_), shape_{ m_,n_ } {}
		constexpr Matrix(Index m_, Index n_, const T& value) requires is_dynamic : data_(m_* n_), shape_{ m_,n_ } {
			fill(value);
		}

		constexpr Matrix(Index m_, Index n_, const std::initializer_list<T> elems) requires (is_dynamic) : Matrix(m_, n_) {
			size_type h = std::min(elems.size(), size());
			std::copy(std::begin(elems), std::begin(elems) + h, begin());
			std::fill(begin() + h, end(), T{});
		}

		constexpr Matrix(Index m_, Index n_, const storage_type& elems) requires is_dynamic : Matrix(m_, n_) {
			std::copy(std::begin(elems), std::end(elems), begin());
		}

		constexpr Matrix(Index m_, Index n_, storage_type&& elems) noexcept requires is_dynamic : shape_{ m_,n_ }, data_{ std::move(elems) } {
			data_.resize(size());
		}


		//
		// Size and access
		// 

		constexpr size_type rows() const noexcept {
			if constexpr (is_dynamic) return shape_.m;
			else return m;
		}

		constexpr size_type cols() const noexcept {
			if constexpr (is_dynamic) return shape_.n;
			else return n;
		}

		constexpr size_type size() const noexcept { return rows() * cols(); }
		constexpr size_type max_size() const noexcept { return size(); }
		constexpr bool empty() const noexcept { return size() == 0; }

		constexpr T* data() noexcept { return data_.data(); }
		constexpr const T* data() const noexcept { return data_.data(); }

		constexpr T& operator()(size_type i, size_type j) noexcept { return data_[j + i * cols()]; }
		constexpr const T& operator()(size_type i, size_type j) const noexcept { return data_[j + i * cols()]; }

		// Bounds checked (raises exception if index is bad)
		constexpr T& at(size_type i, size_type j) { return data_.at(j + i * cols()); }
		constexpr const T& at(size_type i, size_type j) const { return data_.at(j + i * cols()); }

		constexpr void fill(const T& c) { std::fill(begin(), end(), c); }
		constexpr void swap(Matrix& a) noexcept { data_.swap(a.data_); }


		constexpr void resize(Index m_, Index n_) requires (is_dynamic) {
			shape_ = { m_,n_ };
			data_.resize(size());
		}

		// Note: the behaviour is undefined if the m*n != size()
		constexpr void reshape(Index m_, Index n_) requires (is_dynamic) {
			MATRIX_VERIFY(m_ * n_ == size(), "The new shape results in a different size than before", Matrix_shape_error);
			shape_ = { m_,n_ };
		}

		// Note: the behaviour is undefined if the m*n != size()
		constexpr void reshape(Shape<false> new_shape) requires (is_dynamic) {
			MATRIX_VERIFY(new_shape.m * new_shape.n == size(), "The new shape results in a different size than before", Matrix_shape_error);
			shape_ = new_shape;
		}

		constexpr const Shape<!is_dynamic>& shape() const { return shape_; }


		//
		// Arithmetic
		// 

		template<class F> constexpr Matrix& apply(F f) { std::for_each(begin(), end(), f); return *this; }
		template<class F> constexpr Matrix& apply(F f, const T& c) { std::for_each(begin(), end(), [&](T& v) {v = f(v, c); }); return *this; }
		template<class F> constexpr Matrix& apply(F f, const Matrix& v) {
			MATRIX_VERIFY(shape_ == v.shape_, "Cannot operate matrices with non-matching dimensions", Matrix_shape_error);
			std::transform(begin(), end(), v.begin(), begin(), f); return *this;
		}

		constexpr Matrix& operator+=(const T& c) { return apply(std::plus<T>(), c); }
		constexpr Matrix& operator-=(const T& c) { return apply(std::minus<T>(), c); }
		constexpr Matrix& operator*=(const T& c) { return apply(std::multiplies<T>(), c); }
		constexpr Matrix& operator/=(const T& c) { return apply(divides<T>(), c); }
		constexpr Matrix& operator%=(const T& c) { return apply(modulus<T>(), c); }
		constexpr Matrix& operator+=(const Matrix& a) { return apply(std::plus<T>(), a); }
		constexpr Matrix& operator-=(const Matrix& a) { return apply(std::minus<T>(), a); }

		constexpr Matrix operator+(const T& c) const { return Matrix(*this) += c; }
		constexpr Matrix operator-(const T& c) const { return Matrix(*this) -= c; }
		constexpr Matrix operator*(const T& c) const { return Matrix(*this) *= c; }
		constexpr Matrix operator/(const T& c) const { return Matrix(*this) /= c; }
		constexpr Matrix operator%(const T& c) const { return Matrix(*this) %= c; }
		constexpr Matrix operator+(const Matrix& a) const { return Matrix(*this) += a; }
		constexpr Matrix operator-(const Matrix& a) const { return Matrix(*this) -= a; }


		constexpr Matrix<T, n, m> transpose() const {
			Matrix<T, n, m> result{ shape_.swapped() };
			for (size_type i = 0; i < rows(); ++i)
				for (size_type j = 0; j < cols(); ++j)
					result(j, i) = (*this)(i, j);
			return result;
		}

		template<size_type p>
		constexpr Matrix<T, m, p> operator*(const Matrix<T, n, p>& a) const {
			MATRIX_VERIFY(cols() == a.rows(), "Cannot multipliy matrices with non-matching dimensions", Matrix_shape_error);
			Matrix<T, m, p> result(shape() * a.shape());

			for (size_type i = 0; i < rows(); ++i) {
				for (size_type j = 0; j < a.cols(); ++j) {
					T value{};
					for (size_type k = 0; k < a.rows(); ++k)
						value += (*this)(i, k) * a(k, j);
					result(i, j) = value;
				}
			}
			return result;
		}


		template<class OtherScalar, class ResultScalar, size_type p>
		constexpr friend Matrix<T, m, p> operator*(const Matrix<T, m, n>& a, const Matrix<OtherScalar, n, p>& b) {
			MATRIX_VERIFY(a.cols() == b.rows(), "Cannot multipliy matrices with non-matching dimensions", Matrix_shape_error);
			Matrix<ResultScalar, m, p> result(a.shape() * b.shape());

			for (size_type i = 0; i < a.rows(); ++i) {
				for (size_type j = 0; j < b.cols(); ++j) {
					ResultScalar value{};
					for (size_type k = 0; k < a.cols(); ++k)
						value += a(i, k) * b(k, j);
					result(i, j) = value;
				}
			}
			return result;
		}


		//
		// Iterators and block access
		//

		constexpr iterator begin() noexcept { return iterator(data_.data(), 0); }
		constexpr iterator end() noexcept { return iterator(data_.data(), size()); }
		constexpr const_iterator begin() const noexcept { return const_iterator(data_.data(), 0); }
		constexpr const_iterator end() const noexcept { return const_iterator(data_.data(), size()); }
		constexpr const_iterator cbegin() const noexcept { return begin(); }
		constexpr const_iterator cend() const noexcept { return end(); }

		constexpr reverse_iterator rbegin() noexcept { return reverse_iterator(end()); }
		constexpr reverse_iterator rend() noexcept { return reverse_iterator(begin()); }
		constexpr const_reverse_iterator rbegin() const noexcept { return const_reverse_iterator(end()); }
		constexpr const_reverse_iterator rend() const noexcept { return const_reverse_iterator(begin()); }
		constexpr const_reverse_iterator crbegin() const noexcept { return rbegin(); }
		constexpr const_reverse_iterator crend() const noexcept { return rend(); }

		// These are a bit useless with Matrix_view but they could be just a nuance faster. 
		constexpr const_col_iterator col_begin(size_type col) const noexcept { return const_col_iterator(data_.data() + col, cols(), 0); }
		constexpr const_col_iterator col_end(size_type col) const noexcept { return const_col_iterator(data_.data() + col + size(), cols(), 0); }
		constexpr const_row_iterator row_begin(size_type row) const noexcept { return const_row_iterator(data_.data() + row * cols(), 0); }
		constexpr const_row_iterator row_end(size_type row) const noexcept { return const_row_iterator(data_.data() + (row + 1) * cols(), 0); }
		constexpr col_iterator col_begin(size_type col) noexcept { return col_iterator(data_.data() + col, cols(), 0); }
		constexpr col_iterator col_end(size_type col) noexcept { return col_iterator(data_.data() + col + size(), cols(), 0); }
		constexpr row_iterator row_begin(size_type row) noexcept { return row_iterator(data_.data() + row * cols(), 0); }
		constexpr row_iterator row_end(size_type row) noexcept { return row_iterator(data_.data() + (row + 1) * cols(), 0); }


		constexpr Matrix_view<T> block(size_type row, size_type col, size_type rows, size_type cols) {
			MATRIX_VERIFY(row + rows <= this->rows() && col + cols <= this->cols(), "Out of range error at Matrix::block()", Matrix_block_domain_error);
			return Matrix_view<T>{&this->operator()(row, col), this->rows(), this->cols(), rows, cols};
		}
		constexpr Matrix_view<const T> block(size_type row, size_type col, size_type rows, size_type cols) const {
			MATRIX_VERIFY(row + rows <= this->rows() && col + cols <= this->cols(), "Out of range error at Matrix::block()", Matrix_block_domain_error);
			return Matrix_view<const T>{&this->operator()(row, col), this->rows(), this->cols(), rows, cols};
		}
		constexpr Matrix_view<T> row(size_type row) {
			MATRIX_VERIFY(row < rows(), "Out of range error at Matrix::row()", Matrix_block_domain_error);
			return Matrix_view<T>{&this->operator()(row, 0), rows(), cols(), 1, cols()};
		}
		constexpr Matrix_view<const T> row(size_type row) const {
			MATRIX_VERIFY(row < rows(), "Out of range error at Matrix::row()", Matrix_block_domain_error);
			return Matrix_view<const T>{&this->operator()(row, 0), rows(), cols(), 1, cols()};
		}
		constexpr Matrix_view<T> col(size_type col) {
			MATRIX_VERIFY(col < cols(), "Out of range error at Matrix::col()", Matrix_block_domain_error);
			return Matrix_view<T>{&this->operator()(0, col), rows(), cols(), rows(), 1};
		}
		constexpr Matrix_view<const T> col(size_type col) const {
			MATRIX_VERIFY(col < cols(), "Out of range error at Matrix::col()", Matrix_block_domain_error);
			return Matrix_view<const T>{&this->operator()(0, col), rows(), cols(), rows(), 1};
		}


		friend std::ostream& operator<<(std::ostream& os, const Matrix& a) {
			std::vector<size_type> col_widths(a.cols());
			os << std::setfill(' ');
			for (size_type j = 0; j < a.cols(); ++j) {
				size_type max_len = 0;
				for (size_type i = 0; i < a.rows(); ++i) {
					std::ostringstream s;
					s << a(i, j);
					size_type len = s.str().size();
					if (len > max_len) max_len = len;
				}
				col_widths[j] = max_len;
			}

			for (size_type i = 0; i < a.rows(); ++i) {
				os << "| ";
				for (size_type j = 0; j < a.cols(); ++j) {
					std::ostringstream s;
					s << a(i, j);
					os << std::setw(col_widths[j]) << s.str() << ' ';
				}
				os << "|\n";
			}
			os.width(0);
			return os << "\n";
		}


		//
		// Vector operations (specialization for 1×n or m×1 matrices) 
		//
		constexpr bool is_vector() const { return rows() == 1 || cols() == 1; }

		constexpr T& operator[](size_type i) noexcept { return data_[i]; }
		constexpr const T& operator[](size_type i) const noexcept { return data_[i]; }

		constexpr T& x() noexcept requires vector_dimensions_1D<m, n> { return data_[0]; }
		constexpr const T& x() const noexcept requires vector_dimensions_1D<m, n> { return data_[0]; }
		constexpr T& y() noexcept requires vector_dimensions_2D<m, n> { return data_[1]; }
		constexpr const T& y() const noexcept requires vector_dimensions_2D<m, n> { return data_[1]; }
		constexpr T& z() noexcept requires vector_dimensions_3D<m, n> { return data_[2]; }
		constexpr const T& z() const noexcept requires vector_dimensions_3D<m, n> { return data_[2]; }
		constexpr T& w() noexcept requires vector_dimensions_4D<m, n> { return data_[3]; }
		constexpr const T& w() const noexcept requires vector_dimensions_4D<m, n> { return data_[3]; }

		T norm() const requires (vector_dimensions<m, n> || is_dynamic) {
			if constexpr (is_dynamic) {
				MATRIX_VERIFY(is_vector(), "Matrix::norm() is only supported for vectors, not matrices", Matrix_shape_error);
			}
			T inner_product{};
			std::for_each(begin(), end(), [&inner_product](const auto& value) { inner_product += value * value; });
			return std::sqrt(inner_product);
		}

		Matrix& normalize() requires (vector_dimensions<m, n> || is_dynamic) {
			if constexpr (is_dynamic) {
				MATRIX_VERIFY(is_vector(), "Matrix::normalize() is only supported for vectors, not matrices", Matrix_shape_error);
			}
			return *this *= T(1) / norm();
		}

		constexpr T dot(const Matrix& vec) const requires (vector_dimensions<m, n> || is_dynamic) {
			if constexpr (is_dynamic) {
				MATRIX_VERIFY(is_vector() && vec.is_vector(), "Matrix::dot() is only supported for vectors, not matrices", Matrix_shape_error);
			}
			return std::inner_product(begin(), end(), vec.begin(), T{});
		}

		constexpr T operator*(const Matrix& vec) const requires vector_dimensions<m, n>{
			return (*this).dot(vec);
		}

		// Enable cast to T for 1×1 matrix (useful for matrix multiplications that yield a scalar) 
		explicit(false) constexpr operator T() const requires (m == 1 && n == 1) { return data_[0]; }


		//
		// Factory functions
		//

		constexpr static Matrix zero() requires (!is_dynamic) { return Matrix(); }
		constexpr static Matrix zero(Index m_, Index n_) requires (is_dynamic) { return Matrix(m_, n_, 0); }

		constexpr static Matrix identity(Index n_) requires (is_dynamic) {
			Matrix mat(n_, n_);
			for (size_type i = 0; i < n_; ++i) mat(i, i) = T(1);
			return mat;
		}

		constexpr static Matrix identity() requires (m == n && !is_dynamic) {
			Matrix mat;
			for (size_type i = 0; i < n; ++i) mat(i, i) = T(1);
			return mat;
		}

		// Enable cast to other value_type U. T needs to be convertible to U (via static_cast) 
		template<class U> requires std::convertible_to<T, U>
		constexpr Matrix<U, m, n> cast() const {
			Matrix<U, m, n> mat(shape_);
			std::transform(begin(), end(), mat.begin(), [](const T& c) { return static_cast<U>(c); });
			return mat;
		}

		friend constexpr bool operator==(const Matrix& a, const Matrix& b) { return a.data_ == b.data_; }

	private:
		storage_type data_{};
		Shape<!is_dynamic> shape_;

		template<class U> struct divides { constexpr U operator()(const U& l, const U& r) const { return l / r; } };
		template<class U> struct modulus { constexpr U operator()(const U& l, const U& r) const { return l % r; } };
	};


	template<class T>
	Matrix(Index m, Index n, const std::initializer_list<T>& data)->Matrix<T, dynamic, dynamic>;


	template<class T, Index m> using Vector = Matrix<T, m, 1>;
	template<class T, Index n> using RowVector = Matrix<T, 1, n>;



	template<class T, Index m, Index n>
	constexpr bool operator!=(const Matrix<T, m, n>& a, const Matrix<T, m, n>& b) { return !(a == b); }

	template<class T, Index m, Index n>
	constexpr Matrix<T, m, n> operator-(const Matrix<T, m, n>& a) { return a * -1; }

	template<class T, Index m, Index n>
	constexpr Matrix<T, m, n> operator*(const T& c, const Matrix<T, m, n>& a) { return a * c; }

	template<class T, Index m>
	constexpr T distance(const Vector<T, m>& a, const Vector<T, m>& b) { return (a - b).norm(); }

	template<class T, Index m>
	constexpr T distance(const RowVector<T, m>& a, const RowVector<T, m>& b) { return (a - b).norm(); }

	template<class T, Index m, Index n>
	constexpr Matrix<T, m, n> hadamard(const Matrix<T, m, n>& a, const Matrix<T, m, n>& b) { return Matrix<T, m, n>(a).apply(std::multiplies<T>(), b); }


	template<class T, Index n>
	constexpr Matrix<T, n, n> diag(const T(&values)[n]) requires(n != dynamic) {
		Matrix<T, n, n> result;
		for (Index i = 0; i < n; ++i) result(i, i) = values[i];
		return result;
	}

	template<class T, Index n>
	constexpr Matrix<T, n, n> diag(const Vector<T, n>& vec) {
		Matrix<T, n, n> result;
		for (Index i = 0; i < vec.size(); ++i) result(i, i) = vec[i];
		return result;
	}

	template<class T, Index n>
	constexpr Matrix<T, n, n> diag(const RowVector<T, n>& vec) {
		Matrix<T, n, n> result;
		for (Index i = 0; i < n; ++i) result(i, i) = vec[i];
		return result;
	}

	template<class T>
	constexpr Matrix<T, dynamic, dynamic> diag(const Matrix<T>& vec) {
		MATRIX_VERIFY(vec.is_vector(), "Math::diag() needs a vector as argument, not a matrix", Matrix_shape_error);
		const auto n = vec.size();
		Matrix<T> result(n, n);
		for (Index i = 0; i < n; ++i) result(i, i) = vec[i];
		return result;
	}

	template<class T>
	constexpr Matrix<T, dynamic, dynamic> diag(const std::initializer_list<T>& values) {
		Matrix<T> result{ Shape<false>{ values.size(), values.size() } };
		Index i{};
		for (const auto& value : values) {
			result(i, i) = value;
			++i;
		}
		return result;
	}

	template<class T, Index n>
	constexpr Matrix<T, n, n> antidiag(const T(&values)[n]) {
		Matrix<T, n, n> result;
		for (Index i = 0; i < n; ++i) result(n - i - 1, i) = values[i];
		return result;
	}

	template<class T, Index n>
	constexpr Matrix<T, n, n> antidiag(const Vector<T, n>& vec) {
		Matrix<T, n, n> result;
		for (Index i = 0; i < vec.size(); ++i) result(n - i - 1, i) = vec[i];
		return result;
	}

	template<class T, Index n>
	constexpr Matrix<T, n, n> antidiag(const RowVector<T, n>& vec) {
		Matrix<T, n, n> result;
		for (Index i = 0; i < vec.size(); ++i) result(n - i - 1, i) = vec[i];
		return result;
	}

	template<class T>
	constexpr Matrix<T, dynamic, dynamic> antidiag(const Matrix<T>& vec) {
		MATRIX_VERIFY(vec.is_vector(), "Math::antidiag() needs a vector as argument, not a matrix", Matrix_shape_error);
		const auto n = vec.size();
		Matrix<T> result(n, n);
		for (Index i = 0; i < n; ++i) result(n - i - 1, i) = vec[i];
		return result;
	}

	template<class T>
	constexpr Matrix<T, dynamic, dynamic> antidiag(const std::initializer_list<T>& values) {
		const auto n = values.size();
		Matrix<T> result{ Shape<false>{ n, n } };
		Index i{};
		for (const auto& value : values) {
			result(n - i - 1, i) = value;
			++i;
		}
		return result;
	}

	using Matrix2f = Matrix<float, 2, 2>;
	using Matrix3f = Matrix<float, 3, 3>;
	using Matrix4f = Matrix<float, 4, 4>;
	using Vector2f = Vector<float, 2>;
	using Vector3f = Vector<float, 3>;
	using Vector4f = Vector<float, 4>;
	using RowVector2f = RowVector<float, 2>;
	using RowVector3f = RowVector<float, 3>;
	using RowVector4f = RowVector<float, 4>;

	using Matrix2d = Matrix<double, 2, 2>;
	using Matrix3d = Matrix<double, 3, 3>;
	using Matrix4d = Matrix<double, 4, 4>;
	using Vector2d = Vector<double, 2>;
	using Vector3d = Vector<double, 3>;
	using Vector4d = Vector<double, 4>;
	using RowVector2d = RowVector<double, 2>;
	using RowVector3d = RowVector<double, 3>;
	using RowVector4d = RowVector<double, 4>;

}
