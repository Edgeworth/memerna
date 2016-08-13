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
#define USE_HACK_MODEL 1

// Like assert, but can't be disabled.
#define verify_expr(expr, ...) \
  do { \
    if (!(expr)) { \
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
};

void LoadEnergyModelFromDataDir(const std::string& data_dir);
void LoadRandomEnergyModel(energy_t min_energy, energy_t max_energy);
uint32_t EnergyModelChecksum();

std::string sgetline(FILE* fp);
std::string sfmt(const char* fmt, ...);
std::string vsfmt(const char* fmt, va_list l);

}

#endif //MEMERNA_COMMON_H
