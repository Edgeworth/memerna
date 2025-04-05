// Copyright 2025 Eliot Courtney.
#ifndef BACKENDS_COMMON_BASE_BOLTZ_PSEUDOFREE_MODEL_H_
#define BACKENDS_COMMON_BASE_BOLTZ_PSEUDOFREE_MODEL_H_
#include <vector>

#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::base {

class BoltzPseudofreeModel {
 public:
  // Pseudofree energies. Ignored if empty.
  std::vector<BoltzEnergy> paired;
  std::vector<BoltzEnergy> unpaired;
  // Cumulative sum of size N+1 (first element is nothing).
  std::vector<BoltzEnergy> unpaired_cum_log;

  [[nodiscard]] BoltzEnergy Unpaired(int n) const {
    if (unpaired.empty()) return ONE_B;
    return unpaired[n];
  }

  // Inclusive range, unlike pf_unpaired_cum directly.
  [[nodiscard]] BoltzEnergy UnpairedCum(int st, int en) const {
    assert(st <= en + 1);
    if (unpaired.empty()) return ONE_B;
    return exp(unpaired_cum_log[en + 1] - unpaired_cum_log[st]);
  }

  [[nodiscard]] BoltzEnergy Paired(int st, int en) const {
    assert(st <= en + 1);
    if (paired.empty()) return ONE_B;
    return paired[st] * paired[en];
  }

  void Load(std::vector<Energy> pf_paired, std::vector<Energy> pf_unpaired);
  void Verify(const Primary& r) const;
};

}  // namespace mrna::md::base

#endif  // BACKENDS_COMMON_BASE_BOLTZ_PSEUDOFREE_MODEL_H_
