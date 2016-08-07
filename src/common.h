#ifndef MEMERNA_COMMON_H
#define MEMERNA_COMMON_H

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <string>
#include <cstdarg>
#include <cstring>
#include <vector>

// Note that changing these might break tests.
#define COMPUTE_CTDS 0
#define USE_HACK_MODEL 1

// Like assert, but can't be disabled.
#define verify_expr(expr, msg) \
  do { \
    if (!(expr)) { \
      fprintf(stderr, "%s:%d failed: %s\n", __func__, __LINE__, msg); \
      exit(1); \
    } \
  } while(0)

namespace memerna {

typedef int8_t base_t;
typedef std::vector<base_t> rna_t;
typedef int32_t energy_t;

struct folded_rna_t {
  rna_t r;
  std::vector<int> p;
};

// Don't change this value. Plays nice with memset.
const energy_t MAX_E = 0x0F0F0F0F;
const energy_t CAP_E = 0x07070707;
const int HAIRPIN_MIN_SZ = 3;

void Init();

std::string sgetline(FILE* fp);
std::string sfmt(const char* fmt, ...);
std::string vsfmt(const char* fmt, va_list l);

template<typename T, unsigned int K>
struct array3d_t {
  typedef T ArrayType[K];
public:
  array3d_t() : data(nullptr), size(0) {}

  array3d_t(std::size_t size) : data(new T[size * size * K]), size(size) {
    memset(data, MAX_E & 0xFF, sizeof(data[0]) * size * size * K);
  }

  array3d_t(const array3d_t&) = delete;

  array3d_t(array3d_t&& o) {*this = std::move(o);}

  array3d_t& operator=(const array3d_t&) = delete;

  array3d_t& operator=(array3d_t&& o) {
    data = o.data;
    size = o.size;
    o.data = nullptr;
    o.size = 0;
    return *this;
  }

  ~array3d_t() {delete[] data;}

  ArrayType* operator[](std::size_t idx) {
    return reinterpret_cast<ArrayType*>(&data[idx * size * K]);
  }

private:
  T* data;
  std::size_t size;
};

uint32_t EnergyModelChecksum();
}

#endif //MEMERNA_COMMON_H
