// Copyright 2022 Eliot Courtney.
#ifndef API_ENERGY_ENERGY_H_
#define API_ENERGY_ENERGY_H_

#include <memory>
#include <string>

#include "model/ctd.h"
#include "model/energy.h"
#include "util/argparse.h"

namespace mrna {
class Structure;
}  // namespace mrna

namespace mrna::erg {

// NEWMODEL: Add as an option here.
inline Opt BuildOptEnergyModel() {
  return Opt(Opt::ARG)
#if ENERGY_PRECISION == 1
      .Default("t04p1")
      .Choice({"t04p1"});
#elif ENERGY_PRECISION == 2
      .Default("t04p2")
      .Choice({"t04p2", "t04p2full", "t12p2", "t22p2"});
#endif
}

inline const Opt OPT_ENERGY_MODEL =
    BuildOptEnergyModel().LongName("energy-model").ShortName("em").Help("energy model to use");
inline const Opt OPT_MEMERNA_DATA = Opt(Opt::ARG)
                                        .LongName("memerna-data")
                                        .ShortName("md")
                                        .Default("./data/")
                                        .Help("data path for memerna data");
inline const Opt OPT_SEED =
    Opt(Opt::ARG).LongName("seed").Help("seed for random energy model for memerna");
inline const Opt OPT_PAIRED_PSEUDOFREE =
    Opt(Opt::ARG)
        .LongName("pf-paired")
        .Multiple()
        .Help("comma separated energies for paired pseudofree energy");
inline const Opt OPT_UNPAIRED_PSEUDOFREE =
    Opt(Opt::ARG)
        .LongName("pf-unpaired")
        .Multiple()
        .Help("comma separated energies for unpaired pseudofree energy");

void RegisterOptsEnergyModel(ArgParse* args);

std::string ModelPathFromArgParse(const ArgParse& args);
std::string ModelPathFromArgParse(const ArgParse& args, const std::string& model);
std::string ModelPath(const std::string& data_dir, const std::string& model);

struct EnergyResult {
  EnergyResult() = default;
  EnergyResult(Energy energy, Ctds ctd, std::unique_ptr<Structure> struc);
  ~EnergyResult();

  EnergyResult(EnergyResult&&) = default;
  EnergyResult& operator=(EnergyResult&&) = default;
  EnergyResult(const EnergyResult&) = delete;
  EnergyResult& operator=(const EnergyResult&) = delete;

  Energy energy = ZERO_E;
  Ctds ctd;  // May be empty if CTDs were not computed.
  std::unique_ptr<Structure> struc;  // May be nullptr if you didn't ask for a Structure.
};

}  // namespace mrna::erg

#endif  // API_ENERGY_ENERGY_H_
