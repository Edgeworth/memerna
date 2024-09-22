#ifndef BACKENDS_COMMON_BASE_PSEUDOFREE_MODEL_H_
#define BACKENDS_COMMON_BASE_PSEUDOFREE_MODEL_H_

#include <vector>

#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::base {

class PseudofreeModel {
 public:
  // Pseudofree energies. Ignored if empty.
  std::vector<Energy> paired;
  std::vector<Energy> unpaired;
  // Cumulative sum of size N+1 (first element is nothing).
  std::vector<Energy> unpaired_cum;

  [[nodiscard]] constexpr Energy Unpaired(int n) const {
    if (unpaired.empty()) return ZERO_E;
    return unpaired[n];
  }

  // Inclusive range, unlike pf_unpaired_cum directly.
  [[nodiscard]] constexpr Energy UnpairedCum(int st, int en) const {
    if (unpaired.empty()) return ZERO_E;
    return unpaired_cum[en + 1] - unpaired_cum[st];
  }

  [[nodiscard]] constexpr Energy Paired(int st, int en) const {
    if (paired.empty()) return ZERO_E;
    return paired[st] + paired[en];
  }

  void Load(std::vector<Energy> pf_paired, std::vector<Energy> pf_unpaired);
  void Verify(const Primary& r) const;
};

}  // namespace mrna::md::base

#endif  // BACKENDS_COMMON_BASE_PSEUDOFREE_MODEL_H_
