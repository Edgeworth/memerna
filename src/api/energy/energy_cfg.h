// Copyright 2022 Eliot Courtney.
#ifndef API_ENERGY_ENERGY_CFG_H_
#define API_ENERGY_ENERGY_CFG_H_

#include <iosfwd>
#include <string>

#include "util/argparse.h"

namespace mrna::erg {

inline const Opt OPT_SEED =
    Opt(Opt::ARG).LongName("seed").Help("seed for random energy model for memerna");
inline const Opt OPT_MEMERNA_DATA = Opt(Opt::ARG)
                                        .LongName("memerna-data")
                                        .ShortName("md")
                                        .Default("./data/")
                                        .Help("data path for memerna data");

// NEWMODEL: Add as an option here.
inline const Opt OPT_ENERGY_MODEL = Opt(Opt::ARG)
                                        .LongName("energy-model")
                                        .ShortName("em")
#if ENERGY_PRECISION == 1
                                        .Default("t04p1")
                                        .Choice({"t04p1"})
#elif ENERGY_PRECISION == 2
                                        .Default("t04p2")
                                        .Choice({"t04p2", "t04p2full", "t12p2", "t22p2"})
#endif
                                        .Help("energy model to use");
inline const auto OPT_LONELY_PAIRS = mrna::Opt(Opt::ARG)
                                         .LongName("lonely-pairs")
                                         .Choice({"off", "heuristic", "on"})
                                         .Default("heuristic")
                                         .Help("allow lonely pairs");
inline const auto OPT_BULGE_STATES =
    mrna::Opt(Opt::FLAG).LongName("bulge-states").Default(true).Help("count bulge states bonus");
inline const auto OPT_CTD = mrna::Opt(Opt::ARG)
                                .LongName("ctd")
                                .Choice({"none", "d2", "no-coax", "all"})
                                .Default("all")
                                .Help("whether to use CTDs");

void RegisterOpts(ArgParse* args);

std::string ModelPathFromArgParse(const ArgParse& args);
std::string ModelPathFromArgParse(const ArgParse& args, const std::string& model);
std::string ModelPath(const std::string& data_dir, const std::string& model);

struct EnergyCfg {
  enum class LonelyPairs {
    OFF,  // Do not allow lonely pairs.
    HEURISTIC,  // Use a heuristic to disallow lonely pairs (RNAstructure default behaviour).
    ON,  //  Allow lonely pairs.
  };

  enum class Ctd {
    NONE,  //  Do not use CTDs in efn, folding, subopt, partition, etc.
    D2,  // Same as ViennaRNA -d2 in efn, folding, subopt, partition, etc.
    NO_COAX,  //  Use only terminal mismatches and dangling ends in folding, subopt, partition, etc.
    ALL,  //  Use CTDs in folding, subopt, partition, etc.
  };

  // Whether to allow lonely pairs in folding, subopt, partition, etc.
  LonelyPairs lonely_pairs = LonelyPairs::HEURISTIC;

  // Use |bulge_states| to include bonuses for bulge loop states. This is used
  // for minimum free energy like calculations. For partition function like
  // calculations, the states are already handled.
  bool bulge_states = true;

  // TODO(1): Implement and use this.
  // Whether to use CTDs in folding, subopt, partition, etc.
  Ctd ctd = Ctd::ALL;

  static EnergyCfg FromArgParse(const ArgParse& args);
};

std::istream& operator>>(std::istream& str, EnergyCfg::LonelyPairs& o);
std::istream& operator>>(std::istream& str, EnergyCfg::Ctd& o);

}  // namespace mrna::erg

#endif  // API_ENERGY_ENERGY_CFG_H_
