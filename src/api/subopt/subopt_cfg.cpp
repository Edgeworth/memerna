// Copyright 2022 Eliot Courtney.
#include "api/subopt/subopt_cfg.h"

#include <cmath>
#include <ostream>

#include "api/ctx/backend_cfg.h"
#include "util/error.h"

namespace mrna::subopt {

void SuboptCfg::Verify() const {
  verify(strucs >= 0, "subopt strucs must not be negative");
  verify(delta >= ZERO_E, "subopt delta must not be negative");
  verify(delta <= CAP_E, "subopt delta must not exceed CAP");
  verify(std::isfinite(time_secs), "subopt time_secs must be finite");
  verify(time_secs == -1.0F || time_secs >= 0.0F, "subopt time_secs must be -1 or nonnegative");
}

void RegisterOpts(ArgParse* args) {
  RegisterOptsBackendCfg(args);
  args->RegisterOpt(OPT_SUBOPT_DELTA);
  args->RegisterOpt(OPT_SUBOPT_MAX);
  args->RegisterOpt(OPT_SUBOPT_TIME_SECS);
  args->RegisterOpt(OPT_SUBOPT_SORTED);
}

SuboptCfg SuboptCfg::FromArgParse(const ArgParse& args) {
  SuboptCfg cfg;
  args.MaybeSet(OPT_SUBOPT_DELTA, &cfg.delta);
  args.MaybeSet(OPT_SUBOPT_MAX, &cfg.strucs);
  args.MaybeSet(OPT_SUBOPT_TIME_SECS, &cfg.time_secs);
  args.MaybeSet(OPT_SUBOPT_SORTED, &cfg.sorted);
  cfg.Verify();
  return cfg;
}

std::ostream& operator<<(std::ostream& str, const SuboptCfg& o) {
  return str << "SuboptCfg{delta=" << o.delta << ", strucs=" << o.strucs
             << ", time_secs=" << o.time_secs << ", sorted=" << o.sorted << "}";
}

}  // namespace mrna::subopt
