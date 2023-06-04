// Copyright 2023 Eliot Courtney.
#include "api/trace/trace_cfg.h"

#include <ostream>

namespace mrna::trace {

void RegisterOpts(ArgParse* args) { args->RegisterOpt(OPT_TRACE_RANDOM); }

TraceCfg TraceCfg::FromArgParse(const ArgParse& args) {
  TraceCfg cfg;
  args.MaybeSet(OPT_TRACE_RANDOM, &cfg.random);
  return cfg;
}

std::ostream& operator<<(std::ostream& str, const TraceCfg& o) {
  return str << "TraceCfg{"
             << "random=" << o.random << "}";
}

}  // namespace mrna::trace
