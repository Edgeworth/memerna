// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_BASEOPT_ENERGY_BOLTZ_PRECOMP_H_
#define BACKENDS_BASEOPT_ENERGY_BOLTZ_PRECOMP_H_

#include "api/energy/energy_cfg.h"
#include "backends/baseopt/energy/boltz_model.h"
#include "backends/common/base/boltz_precomp_base.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::base::opt {

struct BoltzPrecomp : public BoltzPrecompBase<BoltzModel> {
  erg::EnergyCfg cfg;

  BoltzPrecomp(Primary r, BoltzModel::Ptr bm, erg::EnergyCfg cfg);

  [[nodiscard]] BoltzEnergy Hairpin(int st, int en) const;

  [[nodiscard]] BoltzEnergy TwoLoop(int ost, int oen, int ist, int ien) const;
};

}  // namespace mrna::md::base::opt

#endif  // BACKENDS_BASEOPT_ENERGY_BOLTZ_PRECOMP_H_
