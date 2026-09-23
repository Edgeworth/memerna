// Copyright 2023 Eliot Courtney.
#include "api/trace/trace_cfg.h"

#include <cstdint>
#include <ostream>

#include "util/argparse.h"
#include "util/error.h"

namespace mrna::trace {

void RegisterOpts(ArgParse* args) {
  args->RegisterOpt(OPT_TRACE_RANDOM);
  args->RegisterOpt(OPT_TRACE_SEED);
}

void TraceCfg::Verify() const {
  verify(!seed.has_value() || random, "trace seed requires random traceback");
}

TraceCfg TraceCfg::FromArgParse(const ArgParse& args) {
  TraceCfg cfg;
  args.MaybeSet(OPT_TRACE_RANDOM, &cfg.random);
  cfg.seed = args.MaybeGet<uint_fast32_t>(OPT_TRACE_SEED);
  cfg.Verify();
  return cfg;
}

std::ostream& operator<<(std::ostream& str, const TraceCfg& o) {
  str << "TraceCfg{"
      << "random=" << o.random;
  if (o.seed.has_value()) str << ", seed=" << *o.seed;
  return str << "}";
}

}  // namespace mrna::trace
