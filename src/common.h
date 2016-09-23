#ifndef MEMERNA_COMMON_H
#define MEMERNA_COMMON_H

#include <cassert>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <tuple>
#include <vector>

// Like assert, but can't be disabled.
#define verify_expr(expr, ...)                        \
  do {                                                \
    if (!(expr)) {                                    \
      fprintf(stderr, "%s:%d: ", __func__, __LINE__); \
      fprintf(stderr, __VA_ARGS__);                   \
      fprintf(stderr, "\n");                          \
      exit(1);                                        \
    }                                                 \
  } while (0)

namespace memerna {

typedef int8_t base_t;
typedef std::vector<base_t> primary_t;
typedef int32_t energy_t;

enum Ctd : int8_t {
  CTD_NA,
  CTD_UNUSED,
  CTD_3_DANGLE,
  CTD_5_DANGLE,
  CTD_MISMATCH,
  CTD_LCOAX_WITH_NEXT,
  CTD_LCOAX_WITH_PREV,
  CTD_RCOAX_WITH_NEXT,
  CTD_RCOAX_WITH_PREV,
  CTD_FCOAX_WITH_NEXT,
  CTD_FCOAX_WITH_PREV,
  CTD_SIZE
};

// Don't change this value. Plays nice with memset.
const energy_t MAX_E = 0x0F0F0F0F;
const energy_t CAP_E = 0x07070707;

// -----------------------------------------------
// Values affecting the energy model:
const int HAIRPIN_MIN_SZ = 3;
// N.B. This is for kcal/mol so it's not 8.315.
const double R = 1.985877534e-3;
// This is 37 degrees Celsius. Changing this is not a good idea.
const double T = 310.15;
// Ninio maximum asymmetry.
const energy_t NINIO_MAX_ASYM = 30;
// Maximum size of a twoloop.
const int TWOLOOP_MAX_SZ = 30;

struct secondary_t {
  secondary_t() = default;
  secondary_t(const secondary_t&) = default;
  secondary_t(secondary_t&&) = default;
  secondary_t& operator=(secondary_t&&) = default;

  explicit secondary_t(const primary_t& r_) : r(r_), p(r_.size(), -1) {}
  secondary_t(const primary_t& r_, const std::vector<int>& p_) : r(r_), p(p_) {}

  bool operator==(const secondary_t& o) const { return r == o.r && p == o.p; }
  bool operator!=(const secondary_t& o) const { return !(*this == o); }
  bool operator<(const secondary_t& o) const { return std::tie(r, p) < std::tie(o.r, o.p); }

  primary_t r;
  std::vector<int> p;
};

struct computed_t {
  computed_t() = default;
  computed_t(const computed_t&) = default;
  computed_t(computed_t&&) = default;
  computed_t& operator=(computed_t&&) = default;

  explicit computed_t(const primary_t& r_) : s(r_), base_ctds(r_.size(), CTD_NA), energy(MAX_E) {}

  explicit computed_t(const secondary_t& s_)
      : s(s_), base_ctds(s_.r.size(), CTD_NA), energy(MAX_E) {
    assert(s.r.size() == s.p.size() && s.r.size() == base_ctds.size());
  }

  computed_t(const secondary_t& s_, const std::vector<Ctd>& base_ctds_, energy_t energy_)
      : s(s_), base_ctds(base_ctds_), energy(energy_) {
    assert(s.r.size() == s.p.size() && s.r.size() == base_ctds.size());
  }

  bool operator==(const computed_t& o) const {
    return s == o.s && base_ctds == o.base_ctds && energy == o.energy;
  }
  bool operator!=(const computed_t& o) const { return !(*this == o); }

  secondary_t s;
  std::vector<Ctd> base_ctds;
  energy_t energy;
};

struct computed_energy_comparator_t {
  bool operator()(const computed_t& a, const computed_t& b) const { return a.energy < b.energy; }
};

std::string sgetline(FILE* fp);
std::string sfmt(const char* fmt, ...);
std::string vsfmt(const char* fmt, va_list l);

uint32_t Crc32(const std::string& data);
}

#endif  // MEMERNA_COMMON_H
