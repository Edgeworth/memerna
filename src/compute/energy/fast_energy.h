// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_ENERGY_FAST_ENERGY_H_
#define COMPUTE_ENERGY_FAST_ENERGY_H_

#include <string>
#include <vector>

#include "compute/energy/model.h"
#include "model/globals.h"
#include "model/parsing.h"
#include "util/macros.h"

namespace mrna::energy {

namespace internal {

int MaxNumContiguous(const Primary& r);

}

const int MAX_SPECIAL_HAIRPIN_SZ = 6;

// This is templated because the partition function wants to use it with a different type.
template <typename T, int InitVal>
struct HairpinPrecomp {
  HairpinPrecomp() : num_c(0) {
    std::fill(special, special + sizeof(special) / sizeof(special[0]), InitVal);
  }

  T special[MAX_SPECIAL_HAIRPIN_SZ + 1];
  int num_c;
};

struct Precomp {
  Energy augubranch[4][4];
  Energy min_mismatch_coax;
  Energy min_flush_coax;
  Energy min_twoloop_not_stack;

  std::vector<HairpinPrecomp<Energy, MAX_E>> hairpin;
};

template <typename HairpinPrecomp, typename EM>
std::vector<HairpinPrecomp> PrecomputeHairpin(const Primary& r, const EM& em) {
  std::vector<HairpinPrecomp> pc;
  pc.resize(r.size());
  std::string rna_str = PrimaryToString(r);
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

Precomp PrecomputeData(const Primary& r, const EnergyModel& em);

// Must have global state set.
Energy FastTwoLoop(int ost, int oen, int ist, int ien);
Energy FastHairpin(int st, int en);
inline bool ViableFoldingPair(int st, int en) {
  return CanPair(gr[st], gr[en]) && (en - st - 1 >= HAIRPIN_MIN_SZ) &&
      ((en - st - 3 >= HAIRPIN_MIN_SZ && CanPair(gr[st + 1], gr[en - 1])) ||
          (st > 0 && en < static_cast<int>(gr.size() - 1) && CanPair(gr[st - 1], gr[en + 1])));
}

inline PEnergy Boltzmann(Energy energy) {
  if (energy >= CAP_E) return 0;
  return exp(PEnergy(energy) * (PEnergy(-1) / PEnergy(10.0 * R * T)));
}

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_FAST_ENERGY_H_
