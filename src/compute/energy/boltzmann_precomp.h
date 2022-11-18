// Copyright 2022 Eliot Courtney.
#ifndef COMPUTE_ENERGY_BOLTZMANN_PRECOMP_H_
#define COMPUTE_ENERGY_BOLTZMANN_PRECOMP_H_

#include <memory>
#include <vector>

#include "compute/energy/boltzmann_model.h"
#include "compute/energy/model.h"
#include "compute/energy/precomp.h"
#include "model/constants.h"
#include "model/primary.h"

namespace mrna::energy {

struct BoltzPrecomp {
  BoltzEnergy augubranch[4][4]{};
  std::vector<HairpinPrecomp<BoltzEnergy, -1>> hairpin;

  BoltzPrecomp(Primary r, BoltzEnergyModelPtr bem);

  [[nodiscard]] const EnergyModel& em() const { return bem_->em(); }
  [[nodiscard]] const BoltzEnergyModel& bem() const { return *bem_; }

  [[nodiscard]] BoltzEnergy Hairpin(int st, int en) const;
  [[nodiscard]] BoltzEnergy TwoLoop(int ost, int oen, int ist, int ien) const;

 private:
  Primary r_;
  BoltzEnergyModelPtr bem_;

  void PrecomputeData();
};

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_BOLTZMANN_PRECOMP_H_
