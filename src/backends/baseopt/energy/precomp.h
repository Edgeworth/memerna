// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_BASEOPT_ENERGY_PRECOMP_H_
#define BACKENDS_BASEOPT_ENERGY_PRECOMP_H_

#include "api/energy/energy_cfg.h"
#include "backends/baseopt/energy/model.h"
#include "backends/common/base/precomp_base.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::base::opt {

class Precomp : public PrecompBase<Model> {
 public:
  erg::EnergyCfg cfg;

  Precomp(Primary r, Model::Ptr m, erg::EnergyCfg cfg);

  [[nodiscard]] Energy TwoLoop(int ost, int oen, int ist, int ien) const;

  [[nodiscard]] Energy Hairpin(int st, int en) const;
};

}  // namespace mrna::md::base::opt

#endif  // BACKENDS_BASEOPT_ENERGY_PRECOMP_H_
