// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_COMMON_BASE_BOLTZ_PRECOMP_BASE_H_
#define BACKENDS_COMMON_BASE_BOLTZ_PRECOMP_BASE_H_

#include <utility>
#include <vector>

#include "backends/common/base/precomp_base.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::base {

template <typename BM>
struct BoltzPrecompBase {
  BoltzEnergy augubranch[4][4]{};
  std::vector<HairpinPrecomp<BoltzEnergy>> hairpin;

  BoltzPrecompBase(Primary r, BM::Ptr bm) : r_(std::move(r)), bm_(std::move(bm)) {
    PrecomputeData();
  }

 protected:
  Primary r_;
  BM::Ptr bm_;

 private:
  void PrecomputeData() {
    for (Base i = 0; i < 4; ++i)
      for (Base j = 0; j < 4; ++j) augubranch[i][j] = bm_->multiloop_b * bm_->AuGuPenalty(i, j);
    hairpin = PrecomputeHairpin<HairpinPrecomp<BoltzEnergy>, /*is_boltz=*/true>(r_, *bm_, -1.0);
  }
};

}  // namespace mrna::md::base

#endif  // BACKENDS_COMMON_BASE_BOLTZ_PRECOMP_BASE_H_
