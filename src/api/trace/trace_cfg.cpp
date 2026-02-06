// Copyright 2023 Eliot Courtney.
#include "api/trace/trace_cfg.h"

#include <ostream>

namespace mrna::trace {

void RegisterOpts(ArgParse* args) {
  args->RegisterOpt(OPT_TRACE_RANDOM);
  args->RegisterOpt(OPT_TRACE_SEED);
}

TraceCfg TraceCfg::FromArgParse(const ArgParse& args) {
  TraceCfg cfg;
  args.MaybeSet(OPT_TRACE_RANDOM, &cfg.random);
  cfg.seed = args.MaybeGet<uint_fast32_t>(OPT_TRACE_SEED);
  verify(!cfg.seed.has_value() || cfg.random, "trace-seed requires --trace-random");
  return cfg;
}

std::ostream& operator<<(std::ostream& str, const TraceCfg& o) {
  str << "TraceCfg{"
      << "random=" << o.random;
  if (o.seed.has_value()) str << ", seed=" << *o.seed;
  return str << "}";
}

}  // namespace mrna::trace
