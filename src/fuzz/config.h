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
    Opt().LongName("brute-max").Arg().Help("maximum RNA size to run brute force on");

// MFE fuzzing options:
inline const Opt OPT_FUZZ_MFE = Opt().LongName("mfe").Help("fuzz mfe");
inline const Opt OPT_FUZZ_MFE_RNASTRUCTURE =
    Opt().LongName("mfe-rnastructure").Help("fuzz RNAstructure mfe");
inline const Opt OPT_FUZZ_MFE_BRUTE = Opt().LongName("mfe-brute").Help("fuzz brute mfe");
inline const Opt OPT_FUZZ_MFE_TABLE =
    Opt().LongName("mfe-table").Help("enable checking mfe dp tables");

// Subopt fuzzing:
inline const Opt OPT_FUZZ_SUBOPT = Opt().LongName("subopt").Help("fuzz suboptimal folding");
inline const Opt OPT_FUZZ_SUBOPT_RNASTRUCTURE =
    Opt().LongName("subopt-rnastructure").Help("fuzz RNAstructure suboptimal folding");
inline const Opt OPT_FUZZ_SUBOPT_BRUTE =
    Opt().LongName("subopt-brute").Help("fuzz brute suboptimal folding");
inline const Opt OPT_FUZZ_SUBOPT_MAX =
    Opt()
        .LongName("subopt-max")
        .Arg()
        .Help("maximum number of substructures for subopt max-delta fuzz");
inline const Opt OPT_FUZZ_SUBOPT_DELTA =
    Opt().LongName("subopt-delta").Arg().Help("max energy delta for subopt delta fuzz");

// Partition fuzzing:
inline const Opt OPT_FUZZ_PARTITION = Opt().LongName("partition").Help("fuzz partition function");
inline const Opt OPT_FUZZ_PARTITION_RNASTRUCTURE =
    Opt().LongName("partition-rnastructure").Help("fuzz RNAstructure partition function");
inline const Opt OPT_FUZZ_PARTITION_BRUTE =
    Opt().LongName("partition-brute").Help("fuzz brute partition function");

void RegisterOpts(ArgParse* args);

struct FuzzCfg {
  int brute_max = 30;

  bool mfe = true;
  bool mfe_rnastructure = false;
  bool mfe_brute = true;
  bool mfe_table = true;

  bool subopt = true;
  bool subopt_rnastructure = false;
  bool subopt_brute = true;
  int subopt_max = 5000;
  int subopt_delta = 6;

  bool part = true;
  bool part_rnastructure = false;
  bool part_brute = true;

  std::string Desc();

  static FuzzCfg FromArgParse(const ArgParse& args);
};

}  // namespace mrna::fuzz

#endif  // FUZZ_CONFIG_H_
