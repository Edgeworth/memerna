// Copyright 2021 Eliot Courtney.
#ifndef FUZZ_FUZZ_CFG_H_
#define FUZZ_FUZZ_CFG_H_

#include <string>

#include "api/energy/energy.h"
#include "api/energy/energy_cfg.h"
#include "model/energy.h"
#include "util/argparse.h"

namespace mrna::fuzz {

// Energy model options:
inline const Opt OPT_FUZZ_ENERGY_MODELS =
    erg::BuildOptEnergyModel().LongName("energy-models").ShortName("ems").Multiple();

inline const auto OPT_FUZZ_RANDOM_MODELS =
    mrna::Opt(mrna::Opt::FLAG)
        .LongName("random-models")
        .Help("fuzz random energy models - different each time. c.f. seed");

// Brute force specific options:
// Allows brute force fuzzing to be given a maximum RNA size
inline const Opt OPT_FUZZ_BRUTE_MAX =
    Opt(Opt::ARG).LongName("brute-max").Default(22).Help("maximum RNA size to run brute force on");

// MFE fuzzing options:
inline const Opt OPT_FUZZ_MFE = Opt(Opt::FLAG).LongName("mfe").Default(0).Help("fuzz mfe");
inline const Opt OPT_FUZZ_MFE_RNASTRUCTURE =
    Opt(Opt::FLAG).LongName("mfe-rnastructure").Default(0).Help("fuzz RNAstructure mfe");
inline const Opt OPT_FUZZ_MFE_TABLE =
    Opt(Opt::FLAG).LongName("mfe-table").Default(0).Help("enable checking mfe dp tables");

// Subopt fuzzing:
inline const Opt OPT_FUZZ_SUBOPT =
    Opt(Opt::FLAG).LongName("subopt").Default(0).Help("fuzz suboptimal folding");
inline const Opt OPT_FUZZ_SUBOPT_RNASTRUCTURE = Opt(Opt::FLAG)
                                                    .LongName("subopt-rnastructure")
                                                    .Default(0)
                                                    .Help("fuzz RNAstructure suboptimal folding");
inline const Opt OPT_FUZZ_SUBOPT_STRUCS =
    Opt(Opt::ARG)
        .LongName("subopt-strucs")
        .Default(5000)
        .Help("maximum number of substructures for subopt fuzz");
inline const Opt OPT_FUZZ_SUBOPT_DELTA = Opt(Opt::ARG)
                                             .LongName("subopt-delta")
                                             .Default(E(0.6))
                                             .Help("max energy delta for subopt delta fuzz");

// Partition fuzzing:
inline const Opt OPT_FUZZ_PARTITION =
    Opt(Opt::FLAG).LongName("part").Default(0).Help("fuzz partition function");
inline const Opt OPT_FUZZ_PARTITION_RNASTRUCTURE =
    Opt(Opt::FLAG)
        .LongName("part-rnastructure")
        .Default(0)
        .Help("fuzz RNAstructure partition function");

void RegisterOpts(ArgParse* args);

// Contains all configuration needed to run a fuzzing round.
struct FuzzCfg {
  int brute_max = 22;

  bool mfe = false;
  bool mfe_rnastructure = false;
  bool mfe_table = false;

  bool subopt = false;
  bool subopt_rnastructure = false;
  int subopt_strucs = 5000;
  Energy subopt_delta = E(0.6);

  bool part = false;
  bool part_rnastructure = false;

  // Energy model cfg:
  bool random_models = false;
  erg::EnergyCfg energy_cfg{};
  std::vector<std::string> model_names{};
  std::string data_dir{};

  std::string rnastructure_data_dir{};

  [[nodiscard]] std::string Desc() const;

  static FuzzCfg FromArgParse(const ArgParse& args);
};

}  // namespace mrna::fuzz

#endif  // FUZZ_FUZZ_CFG_H_
