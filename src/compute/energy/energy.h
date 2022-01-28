// Copyright 2022 E.
#ifndef COMPUTE_ENERGY_ENERGY_H_
#define COMPUTE_ENERGY_ENERGY_H_

#include <map>
#include <memory>
#include <string>

#include "model/ctd.h"
#include "model/model.h"
#include "util/argparse.h"

namespace mrna::energy {

inline const Opt OPT_SEED =
    Opt().LongName("seed").Arg().Help("seed for random energy model for memerna");
inline const Opt OPT_MEMERNA_DATA = Opt()
                                        .LongName("memerna-data")
                                        .ShortName("md")
                                        .Arg()
                                        .Help("data path for given energy model for memerna");

void RegisterOpts(ArgParse* args);

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
