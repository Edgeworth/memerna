// Copyright 2016 Eliot Courtney.
#ifndef ENERGY_FAST_ENERGY_H_
#define ENERGY_FAST_ENERGY_H_

#include <string>
#include <vector>

#include "common.h"
#include "energy/energy_model.h"
#include "globals.h"
#include "parsing.h"

namespace mrna {
namespace energy {

namespace internal {

int MaxNumContiguous(const primary_t& r);

}

const int MAX_SPECIAL_HAIRPIN_SZ = 6;

// This is templated because the partition function wants to use it with a different type.
template <typename T, int InitVal>
struct hairpin_precomp_t {
  hairpin_precomp_t() : num_c(0) {
    std::fill(special, special + sizeof(special) / sizeof(special[0]), InitVal);
  }

  T special[MAX_SPECIAL_HAIRPIN_SZ + 1];
  int num_c;
};

struct precomp_t {
  energy_t augubranch[4][4];
  energy_t min_mismatch_coax;
  energy_t min_flush_coax;
  energy_t min_twoloop_not_stack;

  std::vector<hairpin_precomp_t<energy_t, MAX_E>> hairpin;
};

template <typename HairpinPrecomp, typename EM>
std::vector<HairpinPrecomp> PrecomputeHairpin(const primary_t& r, const EM& em) {
  std::vector<HairpinPrecomp> pc;
  pc.resize(r.size());
  std::string rna_str = parsing::PrimaryToString(r);
  for (const auto& hairpinpair : em.hairpin) {
    const auto& str = hairpinpair.first;
    verify(str.size() - 2 <= MAX_SPECIAL_HAIRPIN_SZ, "need to increase MAX_SPECIAL_HAIRPIN_SZ");
    auto pos = rna_str.find(str, 0);
    while (pos != std::string::npos) {
      pc[pos].special[str.size() - 2] = hairpinpair.second;
      pos = rna_str.find(str, pos + 1);
    }
  }
  const int N = static_cast<int>(r.size());
  pc[N - 1].num_c = static_cast<int>(r[N - 1] == C);
  for (int i = N - 2; i >= 0; --i)
    if (r[i] == C) pc[i].num_c = pc[i + 1].num_c + 1;
  return pc;
}

precomp_t PrecomputeData(const primary_t& r, const EnergyModel& em);

// Must have global state set.
energy_t FastTwoLoop(int ost, int oen, int ist, int ien);
energy_t FastHairpin(int st, int en);
inline bool ViableFoldingPair(int st, int en) {
  return CanPair(gr[st], gr[en]) && (en - st - 1 >= HAIRPIN_MIN_SZ) &&
      ((en - st - 3 >= HAIRPIN_MIN_SZ && CanPair(gr[st + 1], gr[en - 1])) ||
          (st > 0 && en < static_cast<int>(gr.size() - 1) && CanPair(gr[st - 1], gr[en + 1])));
}

inline penergy_t Boltzmann(energy_t energy) {
  if (energy >= CAP_E) return 0;
  return exp(penergy_t(energy) * (penergy_t(-1) / penergy_t(10.0 * R * T)));
}

}  // namespace energy
}  // namespace mrna

#endif  // ENERGY_FAST_ENERGY_H_
