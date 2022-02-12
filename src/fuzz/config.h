// Copyright 2021 Eliot Courtney.
#ifndef FUZZ_CONFIG_H_
#define FUZZ_CONFIG_H_

#include <string>

#include "util/argparse.h"

namespace mrna::fuzz {

// Note that fuzz doesn't use flags from anything else - it has all custom flags.

// Brute force specific options:
// Allows brute force fuzzing to be given a maximum RNA size
inline const Opt OPT_FUZZ_BRUTE_MAX =
    Opt(Opt::ARG).LongName("brute-max").Help("maximum RNA size to run brute force on");

// MFE fuzzing options:
inline const Opt OPT_FUZZ_MFE = Opt(Opt::FLAG).LongName("mfe").Help("fuzz mfe");
inline const Opt OPT_FUZZ_MFE_RNASTRUCTURE =
    Opt(Opt::FLAG).LongName("mfe-rnastructure").Help("fuzz RNAstructure mfe");
inline const Opt OPT_FUZZ_MFE_TABLE =
    Opt(Opt::FLAG).LongName("mfe-table").Help("enable checking mfe dp tables");

// Subopt fuzzing:
inline const Opt OPT_FUZZ_SUBOPT =
    Opt(Opt::FLAG).LongName("subopt").Help("fuzz suboptimal folding");
inline const Opt OPT_FUZZ_SUBOPT_RNASTRUCTURE =
    Opt(Opt::FLAG).LongName("subopt-rnastructure").Help("fuzz RNAstructure suboptimal folding");
inline const Opt OPT_FUZZ_SUBOPT_STRUCS =
    Opt(Opt::ARG).LongName("subopt-strucs").Help("maximum number of substructures for subopt fuzz");
inline const Opt OPT_FUZZ_SUBOPT_DELTA =
    Opt(Opt::ARG).LongName("subopt-delta").Help("max energy delta for subopt delta fuzz");

// Partition fuzzing:
inline const Opt OPT_FUZZ_PARTITION =
    Opt(Opt::FLAG).LongName("part").Help("fuzz partition function");
inline const Opt OPT_FUZZ_PARTITION_RNASTRUCTURE =
    Opt(Opt::FLAG).LongName("part-rnastructure").Help("fuzz RNAstructure partition function");

void RegisterOpts(ArgParse* args);

struct FuzzCfg {
  int brute_max = 22;

  bool mfe = true;
  bool mfe_rnastructure = false;
  bool mfe_table = true;

  bool subopt = true;
  bool subopt_rnastructure = false;
  int subopt_strucs = 5000;
  int subopt_delta = 6;

  bool part = true;
  bool part_rnastructure = false;

  [[nodiscard]] std::string Desc() const;

  static FuzzCfg FromArgParse(const ArgParse& args);
};

}  // namespace mrna::fuzz

#endif  // FUZZ_CONFIG_H_
