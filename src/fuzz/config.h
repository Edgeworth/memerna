// Copyright 2021 Eliot Courtney.
#ifndef FUZZ_CONFIG_H_
#define FUZZ_CONFIG_H_

#include <cstdint>
#include <map>
#include <string>

#include "util/argparse.h"

namespace mrna::fuzz {

inline const std::map<std::string, Opt> FUZZ_OPTS = {
    {"random", Opt("use random energy models (disables comparison to RNAstructure)")},
    {"table-check", Opt("enable comparing dp tables between memerna and rnastructure")},
    {"brute-cutoff", Opt("maximum rna size to run brute force on").Arg()},
    {"brute-subopt-max", Opt("maximum number of substructures for brute force fuzz").Arg()},
    {"mfe-rnastructure", Opt("enable rnastructure testing")},
    {"subopt", Opt("enable fuzzing suboptimal folding")},
    {"subopt-rnastructure", Opt("test rnastructure suboptimal folding")},
    {"subopt-max", Opt("maximum number of substructures for subopt max-delta fuzz").Arg()},
    {"subopt-delta", Opt("delta for subopt delta fuzz").Arg()},
    {"partition", Opt("enable fuzzing partition function")},
    {"partition-rnastructure", Opt("test rnastructure partition function")},
};

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
