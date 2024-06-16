// Copyright 2022 Eliot Courtney.
#ifndef API_SUBOPT_SUBOPT_CFG_H_
#define API_SUBOPT_SUBOPT_CFG_H_

#include <limits>

#include "model/energy.h"
#include "util/argparse.h"

namespace mrna::subopt {

inline const auto OPT_SUBOPT_DELTA =
    Opt(Opt::ARG).LongName("subopt-delta").Help("maximum energy delta from minimum");
inline const auto OPT_SUBOPT_MAX =
    Opt(Opt::ARG).LongName("subopt-strucs").Help("maximum number of reported structures");
inline const auto OPT_SUBOPT_TIME_SECS =
    Opt(Opt::ARG)
        .LongName("subopt-time-secs")
        .Help("maximum time in seconds to run suboptimal folding for");
inline const auto OPT_SUBOPT_SORTED = Opt(Opt::FLAG)
                                          .LongName("subopt-sorted")
                                          .Default(true)
                                          .Help("if the structures should be sorted");

void RegisterOpts(ArgParse* args);

struct SuboptCfg {
  constexpr static const int MAX_STRUCTURES = std::numeric_limits<int>::max() / 2;

  Energy delta = CAP_E;  // maximum energy delta from minimum
  int strucs = MAX_STRUCTURES;  // maximum number of structures to report
  float time_secs = -1.0;  // maximum time in seconds to run suboptimal folding for
  bool sorted = true;  // if the structures should be sorted

  static SuboptCfg FromArgParse(const ArgParse& args);
};

}  // namespace mrna::subopt

#endif  // API_SUBOPT_SUBOPT_CFG_H_
