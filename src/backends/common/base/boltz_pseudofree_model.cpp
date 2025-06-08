// Copyright 2025 Eliot Courtney.
#include "backends/common/base/boltz_pseudofree_model.h"

#include <vector>

#include "util/error.h"

namespace mrna::md::base {

void BoltzPseudofreeModel::Load(std::vector<Energy> pf_paired, std::vector<Energy> pf_unpaired) {
  paired.resize(pf_paired.size());
  for (size_t i = 0; i < pf_paired.size(); ++i) paired[i] = pf_paired[i].Boltz();

  unpaired.resize(pf_unpaired.size());
  for (size_t i = 0; i < pf_unpaired.size(); ++i) unpaired[i] = pf_unpaired[i].Boltz();

  if (!unpaired.empty()) {
    unpaired_cum_log.resize(unpaired.size() + 1);
    unpaired_cum_log[0] = ZERO_B;
    for (size_t i = 0; i < unpaired.size(); ++i)
      unpaired_cum_log[i + 1] = unpaired_cum_log[i] + pf_unpaired[i].LogBoltz();
  }
}

void BoltzPseudofreeModel::Verify(const Primary& r) const {
  if (!paired.empty())
    verify(paired.size() == r.size(), "pseudofree paired must be same length as seq");
  if (!unpaired.empty())
    verify(unpaired.size() == r.size(), "pseudofree unpaired must be same length as seq");
}

}  // namespace mrna::md::base
