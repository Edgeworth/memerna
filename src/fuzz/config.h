// Copyright 2021 Eliot Courtney.
#ifndef FUZZ_CONFIG_H_
#define FUZZ_CONFIG_H_

#include <cstdint>

#include "util/argparse.h"

namespace mrna::fuzz {

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
};

inline const std::map<std::string, opt_t> OPTIONS = {
    {"random", opt_t("use random energy models (disables comparison to RNAstructure)")},
    {"table-check", opt_t("enable comparing dp tables between memerna and rnastructure")},
    {"brute-cutoff", opt_t("maximum rna size to run brute force on").Arg()},
    {"brute-subopt-max", opt_t("maximum number of substructures for brute force fuzz").Arg()},
    {"mfe-rnastructure", opt_t("enable rnastructure testing")},
    {"subopt", opt_t("enable fuzzing suboptimal folding")},
    {"subopt-rnastructure", opt_t("test rnastructure suboptimal folding")},
    {"subopt-max", opt_t("maximum number of substructures for subopt max-delta fuzz").Arg()},
    {"subopt-delta", opt_t("delta for subopt delta fuzz").Arg()},
    {"partition", opt_t("enable fuzzing partition function")},
    {"partition-rnastructure", opt_t("test rnastructure partition function")},
};

FuzzCfg FuzzCfgFromArgParse(const ArgParse& args);

}  // namespace mrna::fuzz

#endif  // FUZZ_CONFIG_H_
