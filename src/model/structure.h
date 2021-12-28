// Copyright 2021 Eliot Courtney.
#ifndef MODEL_STRUCTURE_H_
#define MODEL_STRUCTURE_H_

#include <vector>

#include "model/base.h"
#include "util/float.h"
#include "util/macros.h"

namespace mrna {

typedef std::vector<Base> Primary;
typedef int32_t Energy;
typedef flt PEnergy;

// Don't change this value. Plays nice with memset.
inline constexpr Energy MAX_E = 0x0F0F0F0F;
inline constexpr Energy CAP_E = 0x07070707;

// -----------------------------------------------
// Values affecting the energy model:
inline constexpr int HAIRPIN_MIN_SZ = 3;
// N.B. This is for kcal/mol so it's not 8.315.
inline constexpr double R = 1.9872036e-3;
// This is 37 degrees Celsius. Changing this is not a good idea.
inline constexpr double T = 310.15;
// Ninio maximum asymmetry.
inline constexpr Energy NINIO_MAX_ASYM = 30;
// Maximum size of a twoloop.
inline constexpr int TWOLOOP_MAX_SZ = 30;

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

struct Secondary {
  Secondary() = default;
  Secondary(const Secondary&) = default;
  Secondary(Secondary&&) = default;
  Secondary& operator=(Secondary&&) = default;

  explicit Secondary(const Primary& r_) : r(r_), p(r_.size(), -1) {}
  Secondary(const Primary& r_, const std::vector<int>& p_) : r(r_), p(p_) {}

  bool operator==(const Secondary& o) const { return r == o.r && p == o.p; }
  bool operator!=(const Secondary& o) const { return !(*this == o); }
  bool operator<(const Secondary& o) const { return std::tie(r, p) < std::tie(o.r, o.p); }

  Primary r;
  std::vector<int> p;
};

// Secondary structure with MFE information and CTDs.
struct Computed {
  Computed() = default;
  Computed(const Computed&) = default;
  Computed(Computed&&) = default;
  Computed& operator=(Computed&&) = default;

  explicit Computed(const Primary& r_) : s(r_), base_ctds(r_.size(), CTD_NA), energy(MAX_E) {}

  explicit Computed(const Secondary& s_) : s(s_), base_ctds(s_.r.size(), CTD_NA), energy(MAX_E) {
    verify(s.r.size() == s.p.size() && s.r.size() == base_ctds.size(), "bug");
  }

  Computed(const Secondary& s_, const std::vector<Ctd>& base_ctds_, Energy energy_)
      : s(s_), base_ctds(base_ctds_), energy(energy_) {
    verify(s.r.size() == s.p.size() && s.r.size() == base_ctds.size(), "bug");
  }

  bool operator==(const Computed& o) const {
    return s == o.s && base_ctds == o.base_ctds && energy == o.energy;
  }
  bool operator!=(const Computed& o) const { return !(*this == o); }

  Secondary s;
  std::vector<Ctd> base_ctds;
  Energy energy;
};

struct ComputedEnergyCmp {
  bool operator()(const Computed& a, const Computed& b) const { return a.energy < b.energy; }
};

}  // namespace mrna

#endif  // MODEL_STRUCTURE_H_
