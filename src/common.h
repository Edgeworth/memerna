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
typedef std::vector<base_t> primary_t;
typedef int32_t energy_t;

enum Ctd {
  CTD_NA,
  CTD_UNUSED,
  CTD_3_DANGLE,
  CTD_5_DANGLE,
  CTD_TERMINAL_MISMATCH,
  CTD_LEFT_MISMATCH_COAX_WITH_NEXT,
  CTD_LEFT_MISMATCH_COAX_WITH_PREV,
  CTD_RIGHT_MISMATCH_COAX_WITH_NEXT,
  CTD_RIGHT_MISMATCH_COAX_WITH_PREV,
  CTD_FLUSH_COAX_WITH_NEXT,
  CTD_FLUSH_COAX_WITH_PREV,
  CTD_SIZE
};

struct secondary_t {
  secondary_t() = default;
  secondary_t(const secondary_t&) = default;
  secondary_t(secondary_t&&) = default;
  secondary_t& operator=(secondary_t&&) = default;

  explicit secondary_t(const primary_t& r_) : r(r_), p(r_.size(), -1) {}
  secondary_t(const primary_t& r_, const std::vector<int>& p_) : r(r_), p(p_) {}

  primary_t r;
  std::vector<int> p;
};

struct computed_t {
  computed_t() = default;
  computed_t(const computed_t&) = default;
  computed_t(computed_t&&) = default;
  computed_t& operator=(computed_t&&) = default;

  explicit computed_t(const primary_t& r_);
  explicit computed_t(const secondary_t& s_);
  computed_t(const secondary_t& s_, const std::vector<Ctd>& base_ctds_, energy_t energy_);

  secondary_t s;
  std::vector<Ctd> base_ctds;
  energy_t energy;
};

std::string sgetline(FILE* fp);
std::string sfmt(const char* fmt, ...);
std::string vsfmt(const char* fmt, va_list l);

uint32_t Crc32(const std::string& data);

}

#endif //MEMERNA_COMMON_H
