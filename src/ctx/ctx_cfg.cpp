// Copyright 2021 Eliot Courtney.
#include "ctx/ctx_cfg.h"

#include "compute/energy/energy_cfg.h"
#include "compute/subopt/subopt_cfg.h"
#include "util/error.h"

namespace mrna::ctx {

void RegisterOpts(ArgParse* args) {
  energy::RegisterOpts(args);
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
  if (s == "0") {
    o = CtxCfg::DpAlg::ZERO;
  } else if (s == "1") {
    o = CtxCfg::DpAlg::ONE;
  } else if (s == "2") {
    o = CtxCfg::DpAlg::TWO;
  } else if (s == "3") {
    o = CtxCfg::DpAlg::THREE;
  } else if (s == "brute") {
    o = CtxCfg::DpAlg::BRUTE;
  } else {
    error("unknown fold option %s", s.c_str());
  }
  return str;
}

std::istream& operator>>(std::istream& str, CtxCfg::SuboptAlg& o) {
  std::string s;
  str >> s;
  if (s == "0") {
    o = CtxCfg::SuboptAlg::ZERO;
  } else if (s == "1") {
    o = CtxCfg::SuboptAlg::ONE;
  } else if (s == "brute") {
    o = CtxCfg::SuboptAlg::BRUTE;
  } else {
    error("unknown suboptimal option %s", s.c_str());
  }
  return str;
}

std::istream& operator>>(std::istream& str, CtxCfg::PartAlg& o) {
  std::string s;
  str >> s;
  if (s == "0") {
    o = CtxCfg::PartAlg::ZERO;
  } else if (s == "1") {
    o = CtxCfg::PartAlg::ONE;
  } else if (s == "brute") {
    o = CtxCfg::PartAlg::BRUTE;
  } else {
    error("unknown partition option %s", s.c_str());
  }
  return str;
}

}  // namespace mrna::ctx
