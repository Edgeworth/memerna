// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_COMMON_BASE_BOLTZ_PRECOMP_BASE_H_
#define BACKENDS_COMMON_BASE_BOLTZ_PRECOMP_BASE_H_

#include <cstddef>
#include <utility>
#include <vector>

#include "backends/common/base/precomp_base.h"
#include "model/base.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::base {

template <typename BM>
struct BoltzPrecompBase {
  BoltzEnergy augubranch[4][4]{};
  std::vector<HairpinPrecomp<BoltzEnergy>> hairpin;
  std::vector<BoltzEnergy> multiloop_c_pow;

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
    multiloop_c_pow = std::vector<BoltzEnergy>(r_.size() + 1, ONE_B);
    for (std::size_t i = 1; i < multiloop_c_pow.size(); ++i)
      multiloop_c_pow[i] = multiloop_c_pow[i - 1] * bm_->multiloop_c;
    hairpin = PrecomputeHairpin<HairpinPrecomp<BoltzEnergy>, /*is_boltz=*/true>(r_, *bm_, -1.0);
  }
};

}  // namespace mrna::md::base

#endif  // BACKENDS_COMMON_BASE_BOLTZ_PRECOMP_BASE_H_
