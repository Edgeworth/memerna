// Copyright 2022 Eliot Courtney.
#ifndef COMPUTE_BRUTE_BRUTE_CFG_H_
#define COMPUTE_BRUTE_BRUTE_CFG_H_

#include "compute/subopt/subopt_cfg.h"
#include "util/argparse.h"

namespace mrna::brute {

void RegisterOpts(ArgParse* args);

struct BruteCfg {
  bool mfe = false;  // Whether to compute the mfe.
  bool subopt = false;  // Whether to compute suboptimal structures.
  bool part = false;  // Whether to compute the partition function.

  subopt::SuboptCfg subopt_cfg = {};  // Suboptimal structure computation configuration.

  static BruteCfg FromArgParse(const ArgParse& args);
};

}  // namespace mrna::brute

#endif  // COMPUTE_BRUTE_BRUTE_CFG_H_
