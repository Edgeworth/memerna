// Copyright 2016 Eliot Courtney.
#ifndef UTIL_ARRAY_H_
#define UTIL_ARRAY_H_

#include <utility>

#include "model/model.h"

namespace mrna {

template <typename T, unsigned int K>
struct Array3D {
  typedef T ArrayType[K];

 public:
  Array3D() : data_(nullptr), size_(0) {}
  explicit Array3D(std::size_t size, const T init_val = MAX_E)
      : data_(new T[size * size * K]), size_(size) {
    std::fill(data_, data_ + size_ * size_ * K, init_val);
  }
  ~Array3D() { delete[] data_; }

  Array3D(const Array3D&) = delete;
  Array3D& operator=(const Array3D&) = delete;
  Array3D(Array3D&& o) : data_(nullptr), size_(0) { *this = std::move(o); }

  Array3D& operator=(Array3D&& o) {
    delete[] data_;
    data_ = o.data_;
    size_ = o.size_;
    o.data_ = nullptr;
    o.size_ = 0;
    return *this;
  }

  ArrayType* operator[](std::size_t idx) {
    return reinterpret_cast<ArrayType*>(&data_[idx * size_ * K]);
  }

  const ArrayType* operator[](std::size_t idx) const {
    return reinterpret_cast<const ArrayType*>(&data_[idx * size_ * K]);
  }

  std::size_t Size() const { return size_; }

 private:
  T* data_;
  std::size_t size_;
};

template <typename T, unsigned int K>
struct Array2D {
 public:
  Array2D() : data_(nullptr), size_(0) {}
  ~Array2D() { delete[] data_; }

  explicit Array2D(std::size_t size, const T init_val = MAX_E)
      : data_(new T[size * K]), size_(size) {
    std::fill(data_, data_ + size_ * K, init_val);
  }

  Array2D(const Array2D&) = delete;
  Array2D& operator=(const Array2D&) = delete;
  Array2D(Array2D&& o) : data_(nullptr), size_(0) { *this = std::move(o); }

  Array2D& operator=(Array2D&& o) {
    delete[] data_;
    data_ = o.data_;
    size_ = o.size_;
    o.data_ = nullptr;
    o.size_ = 0;
    return *this;
  }

  T* operator[](std::size_t idx) { return &data_[idx * K]; }

  const T* operator[](std::size_t idx) const { return &data_[idx * K]; }

 private:
  T* data_;
  std::size_t size_;
};

}  // namespace mrna

#endif  // UTIL_ARRAY_H_
