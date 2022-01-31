// Copyright 2022 E.
#ifndef COMPUTE_SUBOPT_CONFIG_H_
#define COMPUTE_SUBOPT_CONFIG_H_

#include <limits>
#include <memory>

#include "model/model.h"
#include "util/argparse.h"

namespace mrna::subopt {

inline const auto OPT_SUBOPT_DELTA =
    Opt(Opt::ARG).LongName("subopt-delta").Help("maximum energy delta from minimum");
inline const auto OPT_SUBOPT_MAX =
    Opt(Opt::ARG).LongName("subopt-strucs").Help("maximum number of reported structures");
inline const auto OPT_SUBOPT_SORTED =
    Opt(Opt::FLAG).LongName("subopt-sorted").Help("if the structures should be sorted");

void RegisterOpts(ArgParse* args);

struct SuboptCfg {
  constexpr static const int MAX_STRUCTURES = std::numeric_limits<int>::max() / 2;

  Energy delta = CAP_E;  // maximum energy delta from minimum
  int strucs = MAX_STRUCTURES;  // maximum number of structures to report
  bool sorted = false;  // if the structures should be sorted

  static SuboptCfg FromArgParse(const ArgParse& args);
};

}  // namespace mrna::subopt

#endif  // COMPUTE_SUBOPT_CONFIG_H_
