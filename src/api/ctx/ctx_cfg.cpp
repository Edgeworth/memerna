// Copyright 2021 Eliot Courtney.
#include "api/ctx/ctx_cfg.h"

#include "api/energy/energy.h"
#include "api/subopt/subopt_cfg.h"
#include "api/trace/trace_cfg.h"
#include "util/error.h"

namespace mrna {

void RegisterOpts(ArgParse* args) {
  erg::RegisterOptsEnergyModel(args);
  trace::RegisterOpts(args);
  subopt::RegisterOpts(args);
  args->RegisterOpt(OPT_MFE_ALG);
  args->RegisterOpt(OPT_SUBOPT_ALG);
  args->RegisterOpt(OPT_PART_ALG);
}

CtxCfg CtxCfg::FromArgParse(const ArgParse& args) {
  CtxCfg cfg;
  cfg.mfe_alg = args.Get<MfeAlg>(OPT_MFE_ALG);
  cfg.subopt_alg = args.Get<SuboptAlg>(OPT_SUBOPT_ALG);
  cfg.part_alg = args.Get<PartAlg>(OPT_PART_ALG);
  return cfg;
}

std::istream& operator>>(std::istream& str, CtxCfg::MfeAlg& o) {
  std::string s;
  str >> s;
  if (s == "slowest") {
    o = CtxCfg::MfeAlg::SLOWEST;
  } else if (s == "slow") {
    o = CtxCfg::MfeAlg::SLOW;
  } else if (s == "fastest") {
    o = CtxCfg::MfeAlg::FASTEST;
  } else if (s == "lyngso") {
    o = CtxCfg::MfeAlg::LYNGSO;
  } else if (s == "brute") {
    o = CtxCfg::MfeAlg::BRUTE;
  } else {
    fatal("unknown fold option {}", s);
  }
  return str;
}

std::istream& operator>>(std::istream& str, CtxCfg::SuboptAlg& o) {
  std::string s;
  str >> s;
  if (s == "slowest") {
    o = CtxCfg::SuboptAlg::SLOWEST;
  } else if (s == "fastest") {
    o = CtxCfg::SuboptAlg::FASTEST;
  } else if (s == "brute") {
    o = CtxCfg::SuboptAlg::BRUTE;
  } else {
    fatal("unknown suboptimal option {}", s);
  }
  return str;
}

std::istream& operator>>(std::istream& str, CtxCfg::PartAlg& o) {
  std::string s;
  str >> s;
  if (s == "slowest") {
    o = CtxCfg::PartAlg::SLOWEST;
  } else if (s == "fastest") {
    o = CtxCfg::PartAlg::FASTEST;
  } else if (s == "brute") {
    o = CtxCfg::PartAlg::BRUTE;
  } else {
    fatal("unknown partition option {}", s);
  }
  return str;
}

}  // namespace mrna
