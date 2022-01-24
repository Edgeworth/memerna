// Copyright 2022 Eliot Courtney.
#ifndef COMPUTE_ENERGY_ENERGY_H_
#define COMPUTE_ENERGY_ENERGY_H_

#include "model/ctd.h"
#include "model/model.h"
#include "util/argparse.h"

namespace mrna::energy {

// TODO: this rnastructure specific option should go somewhere else?
inline const std::map<std::string, Opt> ENERGY_OPTS = {
    {"seed", Opt("seed for random energy model for memerna").Arg()},
    {"rnastructure-data", Opt("data path for RNAstructure").Arg()},
    {"memerna-data", Opt("data path for given energy model for memerna").Arg()}};

class Structure;

struct EnergyResult {
  EnergyResult() = default;
  EnergyResult(Energy energy, Ctds ctd, std::unique_ptr<Structure> struc);
  ~EnergyResult();

  EnergyResult(EnergyResult&&) = default;
  EnergyResult& operator=(EnergyResult&&) = default;
  EnergyResult(const EnergyResult&) = delete;
  EnergyResult& operator=(const EnergyResult&) = delete;

  Energy energy = 0;
  Ctds ctd;  // May be empty if CTDs were not computed.
  std::unique_ptr<Structure> struc;  // May be nullptr if you didn't ask for a Structure.
};

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_ENERGY_H_
