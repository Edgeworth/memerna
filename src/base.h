#ifndef MEMERNA_BASE_H
#define MEMERNA_BASE_H

#include <cstdint>
#include <vector>
#include <memory>

#include "common.h"

namespace memerna {

typedef int8_t base_t;
typedef std::vector<base_t> rna_t;
typedef int32_t energy_t;

struct folded_rna_t {
  rna_t r;
  std::vector<int> p;
};

struct array2d_t {
public:
  array2d_t() : data(nullptr), size(0) {}
  array2d_t(std::size_t size) : data(new energy_t[size]), size(size) {}
  array2d_t(const array2d_t&) = delete;
  array2d_t(array2d_t&& o) {*this = std::move(o);}
  array2d_t& operator=(const array2d_t&) = delete;
  array2d_t& operator=(array2d_t&& o) {
    data = o.data;
    size = o.size;
    o.data = nullptr;
    o.size = 0;
    return *this;
  }
  ~array2d_t() {delete[] data;}
  energy_t* operator[](std::size_t idx) {return &data[idx * size];}
private:
  energy_t* data;
  std::size_t size;
};

const base_t A = 0, C = 1, G = 2, U = 3;
const int A_b = 1 << A, C_b = 1 << C, G_b = 1 << G, U_b = 1 << U;

inline bool IsPairOf(base_t a, base_t b, int basesA, int basesB) {
  return (basesA & (1 << a)) && (basesB & 1 << b);
}

inline bool IsUnorderedOf(base_t a, base_t b, int basesA, int basesB) {
  return IsPairOf(a, b, basesA, basesB) || IsPairOf(b, a, basesA, basesB);
}

void Init();

base_t CharToBase(char c);

char BaseToChar(base_t b);

}

#endif //MEMERNA_BASE_H
