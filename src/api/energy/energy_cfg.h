// Copyright 2022 Eliot Courtney.
#ifndef API_ENERGY_ENERGY_CFG_H_
#define API_ENERGY_ENERGY_CFG_H_

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <iosfwd>
#include <string>

#include "model/ctd.h"
#include "util/argparse.h"
#include "util/container.h"
#include "util/string.h"

namespace mrna::erg {

void RegisterOptsEnergyCfg(ArgParse* args);

MAKE_ENUM(EnergyModelKind, T04, T12, T22);

struct EnergyCfg {
  MAKE_NESTED_ENUM(LonelyPairs,
      OFF,  // Do not allow lonely pairs.
      HEURISTIC,  // Use a heuristic to disallow lonely pairs (RNAstructure default behaviour).
      ON  //  Allow lonely pairs.
  );

  MAKE_NESTED_ENUM(Ctd,
      //  Do not use CTDs in efn, folding, subopt, partition, etc. Like "d0".
      NONE,
      // Same as ViennaRNA -d2 in efn, folding, subopt, partition, etc.
      D2,
      //  Use only terminal mismatches and dangling ends in folding, subopt, partition. Like "d1".
      NO_COAX,
      //  Use CTDs in folding, subopt, partition, etc. Like "d3".
      ALL);

  // Whether to allow lonely pairs in folding, subopt, partition, etc.
  LonelyPairs lonely_pairs = LonelyPairs::HEURISTIC;

  // Use `bulge_states` to include bonuses for bulge loop states. This is used
  // for minimum free energy like calculations. For partition function like
  // calculations, the states are already handled.
  bool bulge_states = true;

  // TODO(1): Implement and use this.
  // Whether to use CTDs in folding, subopt, partition, etc.
  Ctd ctd = Ctd::ALL;

  static EnergyCfg FromArgParse(const ArgParse& args);

  [[nodiscard]] constexpr bool UseCoaxialStacking() const { return ctd == Ctd::ALL; }

  [[nodiscard]] constexpr bool UseDangleMismatch() const {
    return ctd == Ctd::ALL || ctd == Ctd::NO_COAX;
  }

  [[nodiscard]] constexpr bool UseD2() const { return ctd == Ctd::D2; }

  [[nodiscard]] auto ParseSeqCtdString(
      const std::string& prim_str, const std::string& ctd_str) const {
    return ::mrna::ParseSeqCtdString(prim_str, ctd_str, UseD2());
  }

  [[nodiscard]] auto ToCtdString(const Secondary& s, const Ctds& ctds) const {
    return ctds.ToString(s, UseD2());
  }
};

// Description of the support of an algorithm for each energy configuration.
struct EnergyCfgSupport {
  smallvec<EnergyCfg::LonelyPairs, 10> lonely_pairs{};
  smallvec<bool, 2> bulge_states{};
  smallvec<EnergyCfg::Ctd, 10> ctd{};

  void VerifySupported(const std::string& name, const EnergyCfg& cfg) const;
};

std::ostream& operator<<(std::ostream& str, const EnergyCfg& o);

inline const auto OPT_LONELY_PAIRS = mrna::Opt(Opt::ARG)
                                         .LongName("lonely-pairs")
                                         .ChoiceEnum<EnergyCfg::LonelyPairs>()
                                         .Default(EnergyCfg::LonelyPairs::HEURISTIC)
                                         .Help("allow lonely pairs");
inline const auto OPT_BULGE_STATES =
    mrna::Opt(Opt::FLAG).LongName("bulge-states").Default(true).Help("count bulge states bonus");
inline const auto OPT_CTD = mrna::Opt(Opt::ARG)
                                .LongName("ctd")
                                .ChoiceEnum<EnergyCfg::Ctd>()
                                .Default(EnergyCfg::Ctd::ALL)
                                .Help("whether to use CTDs");

}  // namespace mrna::erg

template <>
struct fmt::formatter<mrna::erg::EnergyCfg> : ostream_formatter {};

#endif  // API_ENERGY_ENERGY_CFG_H_
