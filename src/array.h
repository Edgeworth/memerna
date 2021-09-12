// Copyright 2016 E.
#ifndef ARRAY_H_
#define ARRAY_H_

#include <utility>

#include "common.h"

namespace mrna {

template <typename T, unsigned int K>
struct array3d_t {
  typedef T ArrayType[K];

 public:
  array3d_t() : data(nullptr), size(0) {}
  explicit array3d_t(std::size_t size_, const T init_val = MAX_E)
      : data(new T[size_ * size_ * K]), size(size_) {
    std::fill(data, data + size * size * K, init_val);
  }
  ~array3d_t() { delete[] data; }

  array3d_t(const array3d_t&) = delete;
  array3d_t& operator=(const array3d_t&) = delete;
  array3d_t(array3d_t&& o) : data(nullptr), size(0) { *this = std::move(o); }

  array3d_t& operator=(array3d_t&& o) {
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
struct array2d_t {
 public:
  array2d_t() : data(nullptr), size(0) {}
  ~array2d_t() { delete[] data; }

  explicit array2d_t(std::size_t size_, const T init_val = MAX_E)
      : data(new T[size_ * K]), size(size_) {
    std::fill(data, data + size * K, init_val);
  }

  array2d_t(const array2d_t&) = delete;
  array2d_t& operator=(const array2d_t&) = delete;
  array2d_t(array2d_t&& o) : data(nullptr), size(0) { *this = std::move(o); }

  array2d_t& operator=(array2d_t&& o) {
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

#endif  // ARRAY_H_
