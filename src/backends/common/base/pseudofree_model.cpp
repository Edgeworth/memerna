// Copyright 2025 Eliot Courtney.
#include "backends/common/base/pseudofree_model.h"

#include <utility>
#include <vector>

#include "util/error.h"

namespace mrna::md::base {

void PseudofreeModel::Load(std::vector<Energy> pf_paired, std::vector<Energy> pf_unpaired) {
  paired = std::move(pf_paired);
  unpaired = std::move(pf_unpaired);
  if (!unpaired.empty()) {
    unpaired_cum.resize(unpaired.size() + 1);
    unpaired_cum[0] = ZERO_E;
    for (size_t i = 0; i < unpaired.size(); ++i)
      unpaired_cum[i + 1] = unpaired_cum[i] + unpaired[i];
  }
}

void PseudofreeModel::Verify(const Primary& r) const {
  if (!paired.empty())
    verify(paired.size() == r.size(), "pseudofree paired must be same length as seq");
  if (!unpaired.empty())
    verify(unpaired.size() == r.size(), "pseudofree unpaired must be same length as seq");
}

}  // namespace mrna::md::base
