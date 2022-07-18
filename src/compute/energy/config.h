// Copyright 2022 Eliot Courtney.
#ifndef COMPUTE_ENERGY_CONFIG_H_
#define COMPUTE_ENERGY_CONFIG_H_

#include <iosfwd>
#include <string>

#include "util/argparse.h"

namespace mrna::energy {

inline const Opt OPT_SEED =
    Opt(Opt::ARG).LongName("seed").Help("seed for random energy model for memerna");
inline const Opt OPT_MEMERNA_DATA = Opt(Opt::ARG)
                                        .LongName("memerna-data")
                                        .ShortName("md")
                                        .Default("./data/model/t04")
                                        .Help("data path for given energy model for memerna");
inline const auto OPT_LONELY_PAIRS =
    mrna::Opt(Opt::FLAG).LongName("lonely-pairs").Default(false).Help("allow lonely pairs");
inline const auto OPT_BULGE_STATES =
    mrna::Opt(Opt::FLAG).LongName("bulge-states").Default(true).Help("count bulge states bonus");
inline const auto OPT_CTD = mrna::Opt(Opt::ARG)
                                .LongName("ctd")
                                .Choice({"none", "d2", "no-coax", "all"})
                                .Default("all")
                                .Help("whether to use CTDs");

void RegisterOpts(ArgParse* args);

struct EnergyCfg {
  enum class Ctd {
    NONE,  //  Do not use CTDs in efn, folding, subopt, partition, etc.
    D2,  // Same as ViennaRNA -d2 in efn, folding, subopt, partition, etc.
    NO_COAX,  //  Use only terminal mismatches and dangling ends in folding, subopt, partition, etc.
    ALL,  //  Use CTDs in folding, subopt, partition, etc.
  };

  // Whether to allow lonely pairs in folding, subopt, partition, etc.
  bool lonely_pairs = false;

  // Use |bulge_states| to include bonuses for bulge loop states. This is used
  // for minimum free energy like calculations. For partition function like
  // calculations, the states are already handled.
  bool bulge_states = true;

  // TODO(1): Implement and use this.
  // Whether to use CTDs in folding, subopt, partition, etc.
  Ctd ctd = Ctd::ALL;

  static EnergyCfg FromArgParse(const ArgParse& args);
};

std::istream& operator>>(std::istream& str, EnergyCfg::Ctd& o);

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_CONFIG_H_
