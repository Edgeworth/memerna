// Copyright 2016 E.
#ifndef UTIL_ARRAY_H_
#define UTIL_ARRAY_H_

#include <utility>

#include "model/structure.h"

namespace mrna {

template <typename T, unsigned int K>
struct Array3D {
  typedef T ArrayType[K];

 public:
  Array3D() : data(nullptr), size(0) {}
  explicit Array3D(std::size_t size_, const T init_val = MAX_E)
      : data(new T[size_ * size_ * K]), size(size_) {
    std::fill(data, data + size * size * K, init_val);
  }
  ~Array3D() { delete[] data; }

  Array3D(const Array3D&) = delete;
  Array3D& operator=(const Array3D&) = delete;
  Array3D(Array3D&& o) : data(nullptr), size(0) { *this = std::move(o); }

  Array3D& operator=(Array3D&& o) {
    delete[] data;
    data = o.data;
    size = o.size;
    o.data = nullptr;
    o.size = 0;
    return *this;
  }

  ArrayType* operator[](std::size_t idx) {
    return reinterpret_cast<ArrayType*>(&data[idx * size * K]);
  }

  const ArrayType* operator[](std::size_t idx) const {
    return reinterpret_cast<const ArrayType*>(&data[idx * size * K]);
  }

  std::size_t Size() const { return size; }

 private:
  T* data;
  std::size_t size;
};

template <typename T, unsigned int K>
struct Array2D {
 public:
  Array2D() : data(nullptr), size(0) {}
  ~Array2D() { delete[] data; }

  explicit Array2D(std::size_t size_, const T init_val = MAX_E)
      : data(new T[size_ * K]), size(size_) {
    std::fill(data, data + size * K, init_val);
  }

  Array2D(const Array2D&) = delete;
  Array2D& operator=(const Array2D&) = delete;
  Array2D(Array2D&& o) : data(nullptr), size(0) { *this = std::move(o); }

  Array2D& operator=(Array2D&& o) {
    delete[] data;
    data = o.data;
    size = o.size;
    o.data = nullptr;
    o.size = 0;
    return *this;
  }

  T* operator[](std::size_t idx) { return &data[idx * K]; }

  const T* operator[](std::size_t idx) const { return &data[idx * K]; }

 private:
  T* data;
  std::size_t size;
};

}  // namespace mrna

#endif  // UTIL_ARRAY_H_
