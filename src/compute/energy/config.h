// Copyright 2022 Eliot Courtney.
#ifndef COMPUTE_ENERGY_CONFIG_H_
#define COMPUTE_ENERGY_CONFIG_H_

#include <iosfwd>
#include <string>

#include "util/argparse.h"

namespace mrna::energy {

inline const Opt OPT_SEED =
    Opt().LongName("seed").Arg().Help("seed for random energy model for memerna");
inline const Opt OPT_MEMERNA_DATA = Opt()
                                        .LongName("memerna-data")
                                        .ShortName("md")
                                        .Arg()
                                        .Help("data path for given energy model for memerna");
inline const auto OPT_LONELY_PAIRS =
    mrna::Opt().LongName("lonely-pairs").Help("allow lonely pairs");
inline const auto OPT_CTD = mrna::Opt()
                                .LongName("ctd")
                                .Choice({"none", "no_coax", "all"})
                                .Default("all")
                                .Help("whether to use CTDs");

void RegisterOpts(ArgParse* args);

// TODO: Implement and use these options.
struct EnergyCfg {
  enum class Ctd {
    NONE,  //  Do not use CTDs in folding, subopt, partition, etc.
    NO_COAX,  //  Use only terminal mismatches and dangling ends in folding, subopt, partition, etc.
    ALL,  //  Use CTDs in folding, subopt, partition, etc.
  };

  bool lonely_pairs = false;  // Whether to allow lonely pairs in folding, subopt, partition, etc.
  Ctd ctd = Ctd::ALL;  // Whether to use CTDs in folding, subopt, partition, etc.

  static EnergyCfg FromArgParse(const ArgParse& args);
};

std::istream& operator>>(std::istream& str, EnergyCfg::Ctd& o);

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_CONFIG_H_
