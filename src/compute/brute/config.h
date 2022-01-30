// Copyright 2022 E.
#ifndef COMPUTE_BRUTE_CONFIG_H_
#define COMPUTE_BRUTE_CONFIG_H_

#include "util/argparse.h"

namespace mrna::brute {

void RegisterOpts(ArgParse* args);

struct BruteCfg {
  bool mfe;  // Whether to compute the mfe.
  bool subopt;  // Whether to compute suboptimal structures.
  bool part;  // Whether to compute the partition function.

  static BruteCfg FromArgParse(const ArgParse& args);
};

}  // namespace mrna::brute

#endif  // COMPUTE_BRUTE_CONFIG_H_
