// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_BASE_ENERGY_BOLTZ_PRECOMP_H_
#define BACKENDS_BASE_ENERGY_BOLTZ_PRECOMP_H_

#include "backends/base/energy/boltz_model.h"
#include "backends/common/base/boltz_precomp_base.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::base {

struct BoltzPrecomp : public BoltzPrecompBase<BoltzModel> {
  BoltzPrecomp(Primary r, BoltzModel::Ptr bm);

  [[nodiscard]] BoltzEnergy Hairpin(int st, int en) const;

  [[nodiscard]] BoltzEnergy TwoLoop(int ost, int oen, int ist, int ien) const;
};

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_ENERGY_BOLTZ_PRECOMP_H_
