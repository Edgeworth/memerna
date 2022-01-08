// Copyright 2022 E.
#ifndef COMPUTE_ENERGY_BOLTZMANN_PRECOMP_H_
#define COMPUTE_ENERGY_BOLTZMANN_PRECOMP_H_

#include "compute/energy/boltzmann_model.h"
#include "compute/energy/model.h"
#include "compute/energy/precomp.h"
#include "model/model.h"

namespace mrna::energy {

struct BoltzPrecomp {
  BoltzEnergy augubranch[4][4];
  std::vector<HairpinPrecomp<BoltzEnergy, -1>> hairpin;

  BoltzPrecomp(Primary r, BoltzEnergyModel bem);

  const BoltzEnergyModel& bem() const { return bem_; }

  BoltzEnergy FastHairpin(int st, int en) const;
  BoltzEnergy FastTwoLoop(int ost, int oen, int ist, int ien) const;

 private:
  const Primary r_;
  const BoltzEnergyModel bem_;

  void PrecomputeData();
};

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_BOLTZMANN_PRECOMP_H_
