// Copyright 2022 Eliot Courtney.
#ifndef API_ENERGY_ENERGY_CFG_H_
#define API_ENERGY_ENERGY_CFG_H_

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <boost/container/small_vector.hpp>
#include <iosfwd>
#include <string>

#include "util/argparse.h"

namespace mrna::erg {

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

void RegisterOptsEnergyCfg(ArgParse* args);

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

// Description of the support of an algorithm for each energy configuration.
struct EnergyCfgSupport {
  boost::container::small_vector<EnergyCfg::LonelyPairs, 10> lonely_pairs{};
  boost::container::small_vector<bool, 2> bulge_states{};
  boost::container::small_vector<EnergyCfg::Ctd, 10> ctd{};

  void VerifySupported(const std::string& name, const EnergyCfg& cfg) const;
};

std::ostream& operator<<(std::ostream& str, const EnergyCfg& o);
std::ostream& operator<<(std::ostream& str, const EnergyCfg::LonelyPairs& o);
std::ostream& operator<<(std::ostream& str, const EnergyCfg::Ctd& o);

std::istream& operator>>(std::istream& str, EnergyCfg::LonelyPairs& o);
std::istream& operator>>(std::istream& str, EnergyCfg::Ctd& o);

}  // namespace mrna::erg

template <>
struct fmt::formatter<mrna::erg::EnergyCfg> : ostream_formatter {};

template <>
struct fmt::formatter<mrna::erg::EnergyCfg::LonelyPairs> : ostream_formatter {};

template <>
struct fmt::formatter<mrna::erg::EnergyCfg::Ctd> : ostream_formatter {};

#endif  // API_ENERGY_ENERGY_CFG_H_
