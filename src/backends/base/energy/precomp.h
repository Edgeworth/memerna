// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_BASE_ENERGY_PRECOMP_H_
#define BACKENDS_BASE_ENERGY_PRECOMP_H_

#include "backends/base/energy/model.h"
#include "backends/common/base/precomp_base.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::base {

class Precomp : public PrecompBase<Model> {
 public:
  Precomp(Primary r, Model::Ptr m);

  [[nodiscard]] Energy TwoLoop(int ost, int oen, int ist, int ien) const;

  [[nodiscard]] Energy Hairpin(int st, int en) const;
};

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_ENERGY_PRECOMP_H_
