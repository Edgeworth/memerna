// Copyright 2022 Eliot Courtney.
#ifndef MODELS_T04_ENERGY_BOLTZ_PRECOMP_H_
#define MODELS_T04_ENERGY_BOLTZ_PRECOMP_H_

#include <memory>
#include <vector>

#include "model/energy.h"
#include "model/primary.h"
#include "models/t04/energy/boltz_model.h"
#include "models/t04/energy/model.h"
#include "models/t04/energy/precomp.h"

namespace mrna::md::t04 {

struct BoltzPrecomp {
  BoltzEnergy augubranch[4][4]{};
  std::vector<HairpinPrecomp<BoltzEnergy>> hairpin;

  BoltzPrecomp(Primary r, BoltzModel::Ptr bem);

  [[nodiscard]] const Model& em() const { return bem_->em(); }
  [[nodiscard]] const BoltzModel& bem() const { return *bem_; }

  [[nodiscard]] BoltzEnergy Hairpin(int st, int en) const;
  [[nodiscard]] BoltzEnergy TwoLoop(int ost, int oen, int ist, int ien) const;

 private:
  Primary r_;
  BoltzModel::Ptr bem_;

  void PrecomputeData();
};

}  // namespace mrna::md::t04

#endif  // MODELS_T04_ENERGY_BOLTZ_PRECOMP_H_
