// Copyright 2022 E.
#ifndef COMPUTE_ENERGY_CONFIG_H_
#define COMPUTE_ENERGY_CONFIG_H_

#include <memory>

#include "util/argparse.h"

namespace mrna::energy {

inline const Opt OPT_SEED =
    Opt().LongName("seed").Arg().Help("seed for random energy model for memerna");
inline const Opt OPT_MEMERNA_DATA = Opt()
                                        .LongName("memerna-data")
                                        .ShortName("md")
                                        .Arg()
                                        .Help("data path for given energy model for memerna");

void RegisterOpts(ArgParse* args);

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_CONFIG_H_
