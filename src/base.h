#ifndef MEMERNA_BASE_H
#define MEMERNA_BASE_H

#include <cstdint>
#include <vector>
#include <memory>
#include <cstring>

#include "common.h"

namespace memerna {

typedef int8_t base_t;
typedef std::vector<base_t> rna_t;
typedef int32_t energy_t;

// Don't change this value. Plays nice with memset.
const energy_t MAX_E = 0x0F0F0F0F;
const energy_t CAP_E = 0x07070707;
const int HAIRPIN_MIN_SZ = 3;

const base_t A = 0, C = 1, G = 2, U = 3;
const int A_b = 1 << A, C_b = 1 << C, G_b = 1 << G, U_b = 1 << U;

struct folded_rna_t {
  rna_t r;
  std::vector<int> p;
};

template<typename T>
struct array2d_t {
public:
  array2d_t() : data(nullptr), size(0) {}
  array2d_t(std::size_t size) : data(new T[size * size]), size(size) {
    memset(data, MAX_E & 0xFF, sizeof(data[0]) * size * size);
  }
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
  T* operator[](std::size_t idx) {return &data[idx * size];}

private:
  T* data;
  std::size_t size;
};


inline bool IsPairOf(base_t a, base_t b, int basesA, int basesB) {
  return (basesA & (1 << a)) && (basesB & 1 << b);
}

inline bool IsUnorderedOf(base_t a, base_t b, int basesA, int basesB) {
  return IsPairOf(a, b, basesA, basesB) || IsPairOf(b, a, basesA, basesB);
}

inline bool CanPair(base_t a, base_t b) {
  int combined = (1 << a) | (1 << b);
  return combined == (G_b | U_b) || combined == (G_b | C_b) || combined == (A_b | U_b);
}

inline bool IsAuGu(base_t a, base_t b) {
  int combined = (1 << a) | (1 << b);
  return combined == (G_b | U_b) || combined == (A_b | U_b);
}

void Init();

base_t CharToBase(char c);

char BaseToChar(base_t b);

}

#endif //MEMERNA_BASE_H
