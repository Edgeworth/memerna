// Copyright 2021 E.
#ifndef FUZZ_CONFIG_H_
#define FUZZ_CONFIG_H_

#include <cstdint>
#include <map>
#include <string>

#include "util/argparse.h"

namespace mrna::fuzz {

inline const Opt OPT_RANDOM =
    Opt().LongName("random").Help("use random energy models (disables comparison to RNAstructure)");
inline const Opt OPT_TABLE_CHECK =
    Opt()
        .LongName("table-check")
        .Help("enable comparing dp tables between memerna and rnastructure");
inline const Opt OPT_BRUTE_CUTOFF =
    Opt().LongName("brute-cutoff").Arg().Help("maximum rna size to run brute force on");
inline const Opt OPT_BRUTE_SUBOPT_MAX =
    Opt()
        .LongName("brute-subopt-max")
        .Arg()
        .Help("maximum number of substructures for brute force fuzz");
inline const Opt OPT_MFE_RNASTRUCTURE =
    Opt().LongName("mfe-rnastructure").Help("enable rnastructure testing");
inline const Opt OPT_SUBOPT = Opt().LongName("subopt").Help("enable fuzzing suboptimal folding");
inline const Opt OPT_SUBOPT_RNASTRUCTURE =
    Opt().LongName("subopt-rnastructure").Help("test rnastructure suboptimal folding");
inline const Opt OPT_SUBOPT_MAX =
    Opt()
        .LongName("subopt-max")
        .Arg()
        .Help("maximum number of substructures for subopt max-delta fuzz");
inline const Opt OPT_SUBOPT_DELTA =
    Opt().LongName("subopt-delta").Arg().Help("delta for subopt delta fuzz");
inline const Opt OPT_PARTITION =
    Opt().LongName("partition").Help("enable fuzzing partition function");
inline const Opt OPT_PARTITION_RNASTRUCTURE =
    Opt().LongName("partition-rnastructure").Help("test rnastructure partition function");

void RegisterOpts(ArgParse* args);

struct FuzzCfg {
  bool random_model = false;
  bool mfe_rnastructure = false;
  bool table_check = true;
  uint_fast32_t seed = 0;

  bool subopt = true;
  bool subopt_rnastructure = false;
  int subopt_max = 5000;
  int subopt_delta = 6;  // Same as RNAstructure default.

  int brute_cutoff = 30;
  int brute_subopt_max = 100000;

  bool partition = true;
  bool partition_rnastructure = false;

  std::string Describe();

  static FuzzCfg FromArgParse(const ArgParse& args);
};

}  // namespace mrna::fuzz

#endif  // FUZZ_CONFIG_H_
