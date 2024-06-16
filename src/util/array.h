// Copyright 2016 Eliot Courtney.
#ifndef UTIL_ARRAY_H_
#define UTIL_ARRAY_H_

#include <algorithm>
#include <utility>

namespace mrna {

// Array types optimised for various memory layouts. The naming convention is
// ArrayXDYS, where there are X dynamically sized dimensions, and Y statically
// sized dimensions. This is useful when storing lots of parallel arrays
// near each other, which is used a lot in RNA DPs. These arrays also assume
// that all dynamic dimensions are the same size, i.e. are square.

template <typename T>
struct Array2D {
 public:
  constexpr Array2D() = default;
  constexpr Array2D(std::size_t size, T init_val) : data_(new T[size * size]), size_(size) {
    std::fill(data_, data_ + size_ * size_, init_val);
  }
  constexpr ~Array2D() { delete[] data_; }

  constexpr Array2D(const Array2D&) = delete;
  constexpr Array2D& operator=(const Array2D&) = delete;
  constexpr Array2D(Array2D&& o) noexcept { *this = std::move(o); }

  constexpr Array2D& operator=(Array2D&& o) noexcept {
    delete[] data_;
    data_ = o.data_;
    size_ = o.size_;
    o.data_ = nullptr;
    o.size_ = 0;
    return *this;
  }

  constexpr T* operator[](std::size_t idx) { return &data_[idx * size_]; }
  constexpr const T* operator[](std::size_t idx) const { return &data_[idx * size_]; }

  [[nodiscard]] constexpr std::size_t size() const { return size_; }
  [[nodiscard]] constexpr bool empty() const { return size_ == 0; }

 private:
  T* data_{nullptr};
  std::size_t size_{0};
};

template <typename T, unsigned int K>
struct Array2D1S {
  using ArrayType = T[K];

 public:
  constexpr Array2D1S() = default;
  constexpr Array2D1S(std::size_t size, T init_val) : data_(new T[size * size * K]), size_(size) {
    std::fill(data_, data_ + size_ * size_ * K, init_val);
  }
  constexpr ~Array2D1S() { delete[] data_; }  // NOLINT - false positive.

  constexpr Array2D1S(const Array2D1S&) = delete;
  constexpr Array2D1S& operator=(const Array2D1S&) = delete;
  constexpr Array2D1S(Array2D1S&& o) noexcept { *this = std::move(o); }

  constexpr Array2D1S& operator=(Array2D1S&& o) noexcept {
    delete[] data_;
    data_ = o.data_;
    size_ = o.size_;
    o.data_ = nullptr;
    o.size_ = 0;
    return *this;
  }

  constexpr ArrayType* operator[](std::size_t idx) {
    return reinterpret_cast<ArrayType*>(&data_[idx * size_ * K]);
  }

  constexpr const ArrayType* operator[](std::size_t idx) const {
    return reinterpret_cast<const ArrayType*>(&data_[idx * size_ * K]);
  }

  [[nodiscard]] constexpr std::size_t size() const { return size_; }
  [[nodiscard]] constexpr bool empty() const { return size_ == 0; }

 private:
  T* data_{nullptr};
  std::size_t size_{0};
};

template <typename T, unsigned int K>
struct Array1D1S {
 public:
  constexpr Array1D1S() = default;
  constexpr ~Array1D1S() { delete[] data_; }

  constexpr Array1D1S(std::size_t size, T init_val) : data_(new T[size * K]), size_(size) {
    std::fill(data_, data_ + size_ * K, init_val);
  }

  constexpr Array1D1S(const Array1D1S&) = delete;
  constexpr Array1D1S& operator=(const Array1D1S&) = delete;
  constexpr Array1D1S(Array1D1S&& o) noexcept { *this = std::move(o); }

  constexpr Array1D1S& operator=(Array1D1S&& o) noexcept {
    delete[] data_;
    data_ = o.data_;
    size_ = o.size_;
    o.data_ = nullptr;
    o.size_ = 0;
    return *this;
  }

  constexpr T* operator[](std::size_t idx) { return &data_[idx * K]; }
  constexpr const T* operator[](std::size_t idx) const { return &data_[idx * K]; }

  [[nodiscard]] constexpr std::size_t size() const { return size_; }
  [[nodiscard]] constexpr bool empty() const { return size_ == 0; }

 private:
  T* data_{nullptr};
  std::size_t size_{0};
};

// View of 2D array inside a 3D array.
template <typename T>
struct Array2DView {
 public:
  constexpr Array2DView() = delete;

  constexpr Array2DView(T* data, std::size_t size) : data_(data), size_(size) {}
  constexpr ~Array2DView() = default;

  // Prevent non-local uses to avoid dangling pointers.
  constexpr Array2DView(const Array2DView&) = delete;
  constexpr Array2DView& operator=(const Array2DView&) = delete;
  constexpr Array2DView(Array2DView&& o) = delete;
  constexpr Array2DView& operator=(Array2DView&& o) = delete;

  constexpr T* operator[](std::size_t idx) { return &data_[idx * size_]; }
  constexpr const T* operator[](std::size_t idx) const { return &data_[idx * size_]; }

  [[nodiscard]] constexpr std::size_t size() const { return size_; }
  [[nodiscard]] constexpr bool empty() const { return size_ == 0; }

 private:
  T* data_{nullptr};
  std::size_t size_{0};
};

template <typename T>
struct Array3D {
 public:
  constexpr Array3D() = default;
  constexpr Array3D(std::size_t size, T init_val) : data_(new T[size * size * size]), size_(size) {
    std::fill(data_, data_ + size_ * size_ * size_, init_val);
  }
  constexpr ~Array3D() { delete[] data_; }

  constexpr Array3D(const Array3D&) = delete;
  constexpr Array3D& operator=(const Array3D&) = delete;
  constexpr Array3D(Array3D&& o) noexcept { *this = std::move(o); }

  constexpr Array3D& operator=(Array3D&& o) noexcept {
    delete[] data_;
    data_ = o.data_;
    size_ = o.size_;
    o.data_ = nullptr;
    o.size_ = 0;
    return *this;
  }

  constexpr Array2DView<T> operator[](std::size_t idx) {
    return Array2DView(data_ + idx * size_ * size_, size_);
  }
  constexpr Array2DView<T> operator[](std::size_t idx) const {
    return Array2DView(data_ + idx * size_ * size_, size_);
  }

  [[nodiscard]] constexpr std::size_t size() const { return size_; }
  [[nodiscard]] constexpr bool empty() const { return size_ == 0; }

 private:
  T* data_{nullptr};
  std::size_t size_{0};
};

}  // namespace mrna

#endif  // UTIL_ARRAY_H_
