// Copyright 2022 E.
#ifndef COMPUTE_SUBOPT_CONFIG_H_
#define COMPUTE_SUBOPT_CONFIG_H_

#include <limits>
#include <memory>

#include "model/model.h"
#include "util/argparse.h"

namespace mrna::subopt {

inline const auto OPT_SUBOPT_DELTA =
    mrna::Opt().LongName("subopt-delta").Arg().Help("maximum energy delta from minimum");
inline const auto OPT_SUBOPT_MAX =
    mrna::Opt().LongName("subopt-strucs").Arg().Help("maximum number of reported structures");
inline const auto OPT_SUBOPT_SORTED =
    mrna::Opt().LongName("subopt-sorted").Help("if the structures should be sorted");

void RegisterOpts(ArgParse* args);

struct SuboptCfg {
  constexpr static const int MAX_STRUCTURES = std::numeric_limits<int>::max() / 2;

  Energy delta = MAX_E;  // maximum energy delta from minimum. -1 means no limit
  int strucs = MAX_STRUCTURES;  // maximum number of structures to report.  -1 means no limit
  bool sorted = false;  // if the structures should be sorted

  static SuboptCfg FromArgParse(const ArgParse& args);
};

}  // namespace mrna::subopt

#endif  // COMPUTE_SUBOPT_CONFIG_H_
