// Copyright 2022 Eliot Courtney.
#ifndef COMPUTE_ENERGY_PRECOMP_H_
#define COMPUTE_ENERGY_PRECOMP_H_
#include <cassert>
#include <cmath>
#include <cstdarg>
#include <cstring>
#include <memory>
#include <string>
#include <unordered_map>

#include "compute/energy/model.h"
#include "model/base.h"
#include "model/model.h"
#include "model/primary.h"
#include "util/argparse.h"

namespace mrna::energy {

namespace internal {

int MaxNumContiguous(const Primary& r);

}

inline constexpr int MAX_SPECIAL_HAIRPIN_SZ = 6;

// TODO: move this?
// This is templated because the partition function wants to use it with a different type.
template <typename T, int InitVal>
struct HairpinPrecomp {
  HairpinPrecomp() : num_c(0) {
    std::fill(special, special + sizeof(special) / sizeof(special[0]), InitVal);
  }

  T special[MAX_SPECIAL_HAIRPIN_SZ + 1];
  int num_c;
};

// TODO: move this?
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

class Precomp {
 public:
  Energy augubranch[4][4];
  Energy min_mismatch_coax;
  Energy min_flush_coax;
  Energy min_twoloop_not_stack;

  std::vector<HairpinPrecomp<Energy, MAX_E>> hairpin;

  Precomp(Primary r, EnergyModel em);

  const EnergyModel& em() const { return em_; }

  Energy TwoLoop(int ost, int oen, int ist, int ien) const;
  Energy Hairpin(int st, int en) const;

 private:
  Primary r_;
  EnergyModel em_;

  void PrecomputeData();
};

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_PRECOMP_H_
