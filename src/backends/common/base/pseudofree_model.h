#ifndef BACKENDS_COMMON_BASE_PSEUDOFREE_MODEL_H_
#define BACKENDS_COMMON_BASE_PSEUDOFREE_MODEL_H_

#include <vector>

#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::base {

class PseudofreeModel {
 public:
  // Pseudofree energies. Ignored if empty.
  std::vector<Energy> pf_paired;
  std::vector<Energy> pf_unpaired;
  // Cumulative sum of size N+1 (first element is nothing).
  std::vector<Energy> pf_unpaired_cum;

  [[nodiscard]] constexpr Energy PfUnpaired(int n) const {
    if (pf_unpaired.empty()) return ZERO_E;
    return pf_unpaired[n];
  }

  // Inclusive range, unlike pf_unpaired_cum directly.
  [[nodiscard]] constexpr Energy PfUnpairedCum(int st, int en) const {
    if (pf_unpaired.empty()) return ZERO_E;
    return pf_unpaired_cum[en + 1] - pf_unpaired_cum[st];
  }

  [[nodiscard]] constexpr Energy PfPaired(int st, int en) const {
    if (pf_paired.empty()) return ZERO_E;
    return pf_paired[st] + pf_paired[en];
  }

  void Load(std::vector<Energy> paired, std::vector<Energy> unpaired);
  void Verify(const Primary& r) const;
};

}  // namespace mrna::md::base

#endif  // BACKENDS_COMMON_BASE_PSEUDOFREE_MODEL_H_
