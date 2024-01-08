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
  if (s == "debug") {
    o = CtxCfg::MfeAlg::DEBUG;
  } else if (s == "opt") {
    o = CtxCfg::MfeAlg::OPT;
  } else if (s == "sparse-opt") {
    o = CtxCfg::MfeAlg::SPARSE_OPT;
  } else if (s == "lyngso-sparse-opt") {
    o = CtxCfg::MfeAlg::LYNGSO_SPARSE_OPT;
  } else if (s == "auto") {
    o = CtxCfg::MfeAlg::AUTO;
  } else if (s == "brute") {
    o = CtxCfg::MfeAlg::BRUTE;
  } else {
    fatal("unknown fold option {}", s);
  }
  return str;
}

std::ostream& operator<<(std::ostream& str, const CtxCfg::MfeAlg& o) {
  switch (o) {
  case CtxCfg::MfeAlg::DEBUG: return str << "debug";
  case CtxCfg::MfeAlg::OPT: return str << "opt";
  case CtxCfg::MfeAlg::SPARSE_OPT: return str << "sparse-opt";
  case CtxCfg::MfeAlg::LYNGSO_SPARSE_OPT: return str << "lyngso-sparse-opt";
  case CtxCfg::MfeAlg::AUTO: return str << "auto";
  case CtxCfg::MfeAlg::BRUTE: return str << "brute";
  }
  unreachable();
}

std::istream& operator>>(std::istream& str, CtxCfg::SuboptAlg& o) {
  std::string s;
  str >> s;
  if (s == "debug") {
    o = CtxCfg::SuboptAlg::DEBUG;
  } else if (s == "iterative") {
    o = CtxCfg::SuboptAlg::ITERATIVE;
  } else if (s == "persistent") {
    o = CtxCfg::SuboptAlg::PERSISTENT;
  } else if (s == "auto") {
    o = CtxCfg::SuboptAlg::AUTO;
  } else if (s == "brute") {
    o = CtxCfg::SuboptAlg::BRUTE;
  } else {
    fatal("unknown suboptimal option {}", s);
  }
  return str;
}

std::ostream& operator<<(std::ostream& str, const CtxCfg::SuboptAlg& o) {
  switch (o) {
  case CtxCfg::SuboptAlg::DEBUG: return str << "debug";
  case CtxCfg::SuboptAlg::ITERATIVE: return str << "iterative";
  case CtxCfg::SuboptAlg::PERSISTENT: return str << "persistent";
  case CtxCfg::SuboptAlg::AUTO: return str << "auto";
  case CtxCfg::SuboptAlg::BRUTE: return str << "brute";
  }
  unreachable();
}

std::istream& operator>>(std::istream& str, CtxCfg::PartAlg& o) {
  std::string s;
  str >> s;
  if (s == "debug") {
    o = CtxCfg::PartAlg::DEBUG;
  } else if (s == "opt") {
    o = CtxCfg::PartAlg::OPT;
  } else if (s == "auto") {
    o = CtxCfg::PartAlg::AUTO;
  } else if (s == "brute") {
    o = CtxCfg::PartAlg::BRUTE;
  } else {
    fatal("unknown partition option {}", s);
  }
  return str;
}

std::ostream& operator<<(std::ostream& str, const CtxCfg::PartAlg& o) {
  switch (o) {
  case CtxCfg::PartAlg::DEBUG: return str << "debug";
  case CtxCfg::PartAlg::OPT: return str << "opt";
  case CtxCfg::PartAlg::AUTO: return str << "auto";
  case CtxCfg::PartAlg::BRUTE: return str << "brute";
  }
  unreachable();
}

}  // namespace mrna
