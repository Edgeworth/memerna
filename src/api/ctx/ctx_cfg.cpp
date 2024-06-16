// Copyright 2021 Eliot Courtney.
#include "api/ctx/ctx_cfg.h"

#include "api/subopt/subopt_cfg.h"
#include "api/trace/trace_cfg.h"

namespace mrna {

void RegisterOpts(ArgParse* args) {
  RegisterOptsBackendCfg(args);
  trace::RegisterOpts(args);
  subopt::RegisterOpts(args);
  args->RegisterOpt(OPT_MFE_ALG);
  args->RegisterOpt(OPT_SUBOPT_ALG);
  args->RegisterOpt(OPT_PFN_ALG);
}

CtxCfg CtxCfg::FromArgParse(const ArgParse& args) {
  CtxCfg cfg;
  cfg.mfe_alg = args.Get<MfeAlg>(OPT_MFE_ALG);
  cfg.subopt_alg = args.Get<SuboptAlg>(OPT_SUBOPT_ALG);
  cfg.pfn_alg = args.Get<PfnAlg>(OPT_PFN_ALG);
  return cfg;
}

}  // namespace mrna
