// Copyright 2022 Eliot Courtney.
#ifndef COMPUTE_ENERGY_T04_BOLTZ_PRECOMP_H_
#define COMPUTE_ENERGY_T04_BOLTZ_PRECOMP_H_

#include <memory>
#include <vector>

#include "compute/energy/t04/boltz_model.h"
#include "compute/energy/t04/precomp.h"
#include "model/primary.h"
#include "compute/energy/t04/model.h"
#include "model/energy.h"

namespace mrna::erg::t04 {

struct BoltzPrecomp {
  BoltzEnergy augubranch[4][4]{};
  std::vector<HairpinPrecomp<BoltzEnergy>> hairpin;

  BoltzPrecomp(Primary r, BoltzModelPtr bem);

  [[nodiscard]] const ModelPtr& em() const { return bem_->em(); }
  [[nodiscard]] const BoltzModel& bem() const { return *bem_; }

  [[nodiscard]] BoltzEnergy Hairpin(int st, int en) const;
  [[nodiscard]] BoltzEnergy TwoLoop(int ost, int oen, int ist, int ien) const;

 private:
  Primary r_;
  BoltzModelPtr bem_;

  void PrecomputeData();
};

}  // namespace mrna::erg::t04

#endif  // COMPUTE_ENERGY_T04_BOLTZ_PRECOMP_H_
