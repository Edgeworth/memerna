// Copyright 2022 Eliot Courtney.
#ifndef API_BRUTE_BRUTE_CFG_H_
#define API_BRUTE_BRUTE_CFG_H_

#include "api/subopt/subopt_cfg.h"
#include "util/argparse.h"

namespace mrna::brute {

void RegisterOpts(ArgParse* args);

struct BruteCfg {
  bool mfe = false;  // Whether to compute the mfe.
  bool subopt = false;  // Whether to compute suboptimal structures.
  bool pfn = false;  // Whether to compute the partition function.

  subopt::SuboptCfg subopt_cfg = {};  // Subopt structure computation configuration.

  static BruteCfg FromArgParse(const ArgParse& args);
};

}  // namespace mrna::brute

#endif  // API_BRUTE_BRUTE_CFG_H_
