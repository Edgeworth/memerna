// Copyright 2021 Eliot Courtney.
#include "api/ctx/ctx_cfg.h"

#include "api/energy/energy_cfg.h"
#include "api/subopt/subopt_cfg.h"
#include "util/error.h"

namespace mrna {

void RegisterOpts(ArgParse* args) {
  erg::RegisterOpts(args);
  subopt::RegisterOpts(args);
  args->RegisterOpt(OPT_DP_ALG);
  args->RegisterOpt(OPT_SUBOPT_ALG);
  args->RegisterOpt(OPT_PART_ALG);
}

CtxCfg CtxCfg::FromArgParse(const ArgParse& args) {
  CtxCfg cfg;
  cfg.dp_alg = args.Get<DpAlg>(OPT_DP_ALG);
  cfg.subopt_alg = args.Get<SuboptAlg>(OPT_SUBOPT_ALG);
  cfg.part_alg = args.Get<PartAlg>(OPT_PART_ALG);
  return cfg;
}

std::istream& operator>>(std::istream& str, CtxCfg::DpAlg& o) {
  std::string s;
  str >> s;
  if (s == "slowest") {
    o = CtxCfg::DpAlg::SLOWEST;
  } else if (s == "slow") {
    o = CtxCfg::DpAlg::SLOW;
  } else if (s == "fastest") {
    o = CtxCfg::DpAlg::FASTEST;
  } else if (s == "lyngso") {
    o = CtxCfg::DpAlg::LYNGSO;
  } else if (s == "brute") {
    o = CtxCfg::DpAlg::BRUTE;
  } else {
    error("unknown fold option {}", s);
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
    error("unknown suboptimal option {}", s);
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
    error("unknown partition option {}", s);
  }
  return str;
}

}  // namespace mrna
