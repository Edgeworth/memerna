// Copyright 2023 Eliot Courtney.;
#ifndef API_TRACE_TRACE_CFG_H_
#define API_TRACE_TRACE_CFG_H_

#include <fmt/ostream.h>

#include <iosfwd>
#include <string>

#include "util/argparse.h"

namespace mrna::trace {

inline const auto OPT_TRACE_RANDOM = mrna::Opt(Opt::FLAG)
                                         .LongName("trace-random")
                                         .Default(false)
                                         .Help("take a random MFE trace instead of arbitrary");

void RegisterOpts(ArgParse* args);

struct TraceCfg {
  bool random = false;

  static TraceCfg FromArgParse(const ArgParse& args);
};

std::ostream& operator<<(std::ostream& str, const TraceCfg& o);

}  // namespace mrna::trace

template <>
struct fmt::formatter<mrna::trace::TraceCfg> : ostream_formatter {};

#endif  // API_TRACE_TRACE_CFG_H_
