#include "backends/common/base/pseudofree_model.h"

#include "util/error.h"

namespace mrna::md::base {

void PseudofreeModel::Load(std::vector<Energy> paired, std::vector<Energy> unpaired) {
  pf_paired = std::move(paired);
  pf_unpaired = std::move(unpaired);
  if (!pf_unpaired.empty()) {
    pf_unpaired_cum.resize(pf_unpaired.size() + 1);
    pf_unpaired_cum[0] = ZERO_E;
    for (size_t i = 0; i < pf_unpaired.size(); ++i)
      pf_unpaired_cum[i + 1] = pf_unpaired_cum[i] + pf_unpaired[i];
  }
}

void PseudofreeModel::Verify(const Primary& r) const {
  if (!pf_paired.empty())
    verify(pf_paired.size() == r.size(), "pseudofree paired must be same length as seq");
  if (!pf_unpaired.empty())
    verify(pf_unpaired.size() == r.size(), "pseudofree unpaired must be same length as seq");
}

}  // namespace mrna::md::base
