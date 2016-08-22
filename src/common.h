#ifndef MEMERNA_COMMON_H
#define MEMERNA_COMMON_H

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <string>
#include <cstdarg>
#include <cstring>
#include <vector>

// Like assert, but can't be disabled.
#define verify_expr(expr, ...) \
  do { \
    if (!(expr)) { \
      fprintf(stderr, "%s:%d: ", __func__, __LINE__); \
      fprintf(stderr, __VA_ARGS__); \
      fprintf(stderr, "\n"); \
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
  energy_t energy;
};

std::string sgetline(FILE* fp);
std::string sfmt(const char* fmt, ...);
std::string vsfmt(const char* fmt, va_list l);

uint32_t Crc32(const std::string& data);

}

#endif //MEMERNA_COMMON_H
