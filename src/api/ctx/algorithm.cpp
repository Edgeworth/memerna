// Copyright 2024 Eliot Courtney.
#include "api/ctx/algorithm.h"

#include <fmt/core.h>
#include <iterator>

#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "api/subopt/subopt_cfg.h"
#include "backends/base/energy/model.h"
#include "backends/base/mfe/mfe_debug.h"
#include "backends/base/mfe/mfe_lyngso_sparse_opt.h"
#include "backends/base/mfe/mfe_opt.h"
#include "backends/base/mfe/mfe_sparse_opt.h"
#include "backends/base/pfn/pfn_debug.h"
#include "backends/base/pfn/pfn_opt.h"
#include "backends/base/subopt/subopt_debug.h"
#include "backends/base/subopt/subopt_iterative.h"
#include "backends/base/subopt/subopt_persistent.h"
#include "backends/baseopt/energy/model.h"
#include "backends/baseopt/mfe/mfe_debug.h"
#include "backends/baseopt/mfe/mfe_lyngso_sparse_opt.h"
#include "backends/baseopt/mfe/mfe_opt.h"
#include "backends/baseopt/mfe/mfe_sparse_opt.h"
#include "backends/baseopt/pfn/pfn_debug.h"
#include "backends/baseopt/pfn/pfn_opt.h"
#include "backends/baseopt/subopt/subopt_debug.h"
#include "backends/baseopt/subopt/subopt_iterative.h"
#include "backends/baseopt/subopt/subopt_persistent.h"
#include "backends/stack/energy/model.h"
#include "backends/stack/mfe/mfe_opt.h"
#include "backends/stack/subopt/subopt_iterative.h"
#include "backends/stack/subopt/subopt_persistent.h"
#include "backends/brute/brute.h"
#include "util/error.h"

namespace mrna {

namespace {

constexpr MfeAlg MFE_PRIORITY_BASE[] = {
    MfeAlg::SPARSE_OPT,
    MfeAlg::LYNGSO_SPARSE_OPT,
    MfeAlg::OPT,
    MfeAlg::DEBUG,
};

constexpr MfeAlg MFE_PRIORITY_BASEOPT[] = {
    MfeAlg::SPARSE_OPT,
    MfeAlg::LYNGSO_SPARSE_OPT,
    MfeAlg::OPT,
    MfeAlg::DEBUG,
};

constexpr MfeAlg MFE_PRIORITY_STACK[] = {
    MfeAlg::OPT,
};

constexpr SuboptAlg SUBOPT_PRIORITY_BASE[] = {
    SuboptAlg::ITERATIVE,
    SuboptAlg::PERSISTENT,
    SuboptAlg::ITERATIVE_LOWMEM,
    SuboptAlg::PERSISTENT_LOWMEM,
    SuboptAlg::DEBUG,
};

constexpr SuboptAlg SUBOPT_PRIORITY_BASEOPT[] = {
    SuboptAlg::ITERATIVE,
    SuboptAlg::PERSISTENT,
    SuboptAlg::ITERATIVE_LOWMEM,
    SuboptAlg::PERSISTENT_LOWMEM,
    SuboptAlg::DEBUG,
};

constexpr SuboptAlg SUBOPT_PRIORITY_STACK[] = {
    SuboptAlg::ITERATIVE,
    SuboptAlg::PERSISTENT,
    SuboptAlg::ITERATIVE_LOWMEM,
    SuboptAlg::PERSISTENT_LOWMEM,
};

constexpr PfnAlg PFN_PRIORITY_BASE[] = {
    PfnAlg::OPT,
    PfnAlg::DEBUG,
};

constexpr PfnAlg PFN_PRIORITY_BASEOPT[] = {
    PfnAlg::OPT,
    PfnAlg::DEBUG,
};

}  // namespace

bool BackendIsSupported(BackendKind kind, const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf,
    std::string* reason) {
  switch (kind) {
  case BackendKind::BASE: return md::base::Model::IsSupported(cfg, pf, reason);
  case BackendKind::BASEOPT: return md::base::opt::Model::IsSupported(cfg, pf, reason);
  case BackendKind::STACK: return md::stack::Model::IsSupported(cfg, pf, reason);
  }
  unreachable();
}

smallvec<MfeAlg, EnumCount<MfeAlg>()> MfePriorityForBackend(BackendKind kind, bool include_brute) {
  smallvec<MfeAlg, EnumCount<MfeAlg>()> result;
  switch (kind) {
  case BackendKind::BASE: result.assign(std::begin(MFE_PRIORITY_BASE), std::end(MFE_PRIORITY_BASE)); break;
  case BackendKind::BASEOPT: result.assign(std::begin(MFE_PRIORITY_BASEOPT), std::end(MFE_PRIORITY_BASEOPT)); break;
  case BackendKind::STACK: result.assign(std::begin(MFE_PRIORITY_STACK), std::end(MFE_PRIORITY_STACK)); break;
  }
  if (include_brute) result.push_back(MfeAlg::BRUTE);
  return result;
}

smallvec<SuboptAlg, EnumCount<SuboptAlg>()> SuboptPriorityForBackend(BackendKind kind, bool include_brute) {
  smallvec<SuboptAlg, EnumCount<SuboptAlg>()> result;
  switch (kind) {
  case BackendKind::BASE: result.assign(std::begin(SUBOPT_PRIORITY_BASE), std::end(SUBOPT_PRIORITY_BASE)); break;
  case BackendKind::BASEOPT: result.assign(std::begin(SUBOPT_PRIORITY_BASEOPT), std::end(SUBOPT_PRIORITY_BASEOPT)); break;
  case BackendKind::STACK: result.assign(std::begin(SUBOPT_PRIORITY_STACK), std::end(SUBOPT_PRIORITY_STACK)); break;
  }
  if (include_brute) result.push_back(SuboptAlg::BRUTE);
  return result;
}

smallvec<PfnAlg, EnumCount<PfnAlg>()> PfnPriorityForBackend(BackendKind kind, bool include_brute) {
  smallvec<PfnAlg, EnumCount<PfnAlg>()> result;
  switch (kind) {
  case BackendKind::BASE: result.assign(std::begin(PFN_PRIORITY_BASE), std::end(PFN_PRIORITY_BASE)); break;
  case BackendKind::BASEOPT: result.assign(std::begin(PFN_PRIORITY_BASEOPT), std::end(PFN_PRIORITY_BASEOPT)); break;
  case BackendKind::STACK: break;
  }
  if (include_brute) result.push_back(PfnAlg::BRUTE);
  return result;
}

bool MfeAlgIsSupported(BackendKind kind, MfeAlg alg, const erg::EnergyCfg& cfg,
    const erg::PseudofreeCfg& pf, std::string* reason) {
  if (alg == MfeAlg::AUTO) {
    if (reason) *reason = "AUTO is not a concrete algorithm";
    return false;
  }
  if (!BackendIsSupported(kind, cfg, pf, reason)) return false;
  if (alg == MfeAlg::BRUTE) return md::brute::Brute::IsSupported(cfg, pf, reason);

  switch (kind) {
  case BackendKind::BASE:
    switch (alg) {
    case MfeAlg::AUTO:
    case MfeAlg::BRUTE: break;
    case MfeAlg::DEBUG: return md::base::MfeDebug::IsSupported(cfg, pf, reason);
    case MfeAlg::OPT: return md::base::MfeOpt::IsSupported(cfg, pf, reason);
    case MfeAlg::SPARSE_OPT: return md::base::MfeSparseOpt::IsSupported(cfg, pf, reason);
    case MfeAlg::LYNGSO_SPARSE_OPT:
      return md::base::MfeLyngsoSparseOpt::IsSupported(cfg, pf, reason);
    }
    break;
  case BackendKind::BASEOPT:
    switch (alg) {
    case MfeAlg::AUTO:
    case MfeAlg::BRUTE: break;
    case MfeAlg::DEBUG: return md::base::opt::MfeDebug::IsSupported(cfg, pf, reason);
    case MfeAlg::OPT: return md::base::opt::MfeOpt::IsSupported(cfg, pf, reason);
    case MfeAlg::SPARSE_OPT: return md::base::opt::MfeSparseOpt::IsSupported(cfg, pf, reason);
    case MfeAlg::LYNGSO_SPARSE_OPT:
      return md::base::opt::MfeLyngsoSparseOpt::IsSupported(cfg, pf, reason);
    }
    break;
  case BackendKind::STACK:
    switch (alg) {
    case MfeAlg::AUTO:
    case MfeAlg::BRUTE:
    case MfeAlg::DEBUG:
    case MfeAlg::SPARSE_OPT:
    case MfeAlg::LYNGSO_SPARSE_OPT: break;
    case MfeAlg::OPT: return md::stack::MfeOpt::IsSupported(cfg, pf, reason);
    }
    break;
  }
  if (reason) *reason = fmt::format("not available for {} backend", kind);
  return false;
}

bool SuboptAlgIsSupported(BackendKind kind, SuboptAlg alg, const erg::EnergyCfg& cfg,
    const erg::PseudofreeCfg& pf, const subopt::SuboptCfg& subopt_cfg, std::string* reason) {
  if (alg == SuboptAlg::AUTO) {
    if (reason) *reason = "AUTO is not a concrete algorithm";
    return false;
  }
  if (!BackendIsSupported(kind, cfg, pf, reason)) return false;
  if (alg == SuboptAlg::BRUTE) return md::brute::Brute::IsSupported(cfg, pf, reason);

  switch (kind) {
  case BackendKind::BASE:
    switch (alg) {
    case SuboptAlg::AUTO:
    case SuboptAlg::BRUTE: break;
    case SuboptAlg::DEBUG: return md::base::SuboptDebug::IsSupported(cfg, pf, subopt_cfg, reason);
    case SuboptAlg::ITERATIVE:
      return md::base::SuboptIterative<false>::IsSupported(cfg, pf, subopt_cfg, reason);
    case SuboptAlg::ITERATIVE_LOWMEM:
      return md::base::SuboptIterative<true>::IsSupported(cfg, pf, subopt_cfg, reason);
    case SuboptAlg::PERSISTENT:
      return md::base::SuboptPersistent<false>::IsSupported(cfg, pf, subopt_cfg, reason);
    case SuboptAlg::PERSISTENT_LOWMEM:
      return md::base::SuboptPersistent<true>::IsSupported(cfg, pf, subopt_cfg, reason);
    }
    break;
  case BackendKind::BASEOPT:
    switch (alg) {
    case SuboptAlg::AUTO:
    case SuboptAlg::BRUTE: break;
    case SuboptAlg::DEBUG:
      return md::base::opt::SuboptDebug::IsSupported(cfg, pf, subopt_cfg, reason);
    case SuboptAlg::ITERATIVE:
      return md::base::opt::SuboptIterative<false>::IsSupported(cfg, pf, subopt_cfg, reason);
    case SuboptAlg::ITERATIVE_LOWMEM:
      return md::base::opt::SuboptIterative<true>::IsSupported(cfg, pf, subopt_cfg, reason);
    case SuboptAlg::PERSISTENT:
      return md::base::opt::SuboptPersistent<false>::IsSupported(cfg, pf, subopt_cfg, reason);
    case SuboptAlg::PERSISTENT_LOWMEM:
      return md::base::opt::SuboptPersistent<true>::IsSupported(cfg, pf, subopt_cfg, reason);
    }
    break;
  case BackendKind::STACK:
    switch (alg) {
    case SuboptAlg::AUTO:
    case SuboptAlg::BRUTE:
    case SuboptAlg::DEBUG: break;
    case SuboptAlg::ITERATIVE:
      return md::stack::SuboptIterative<false>::IsSupported(cfg, pf, subopt_cfg, reason);
    case SuboptAlg::ITERATIVE_LOWMEM:
      return md::stack::SuboptIterative<true>::IsSupported(cfg, pf, subopt_cfg, reason);
    case SuboptAlg::PERSISTENT:
      return md::stack::SuboptPersistent<false>::IsSupported(cfg, pf, subopt_cfg, reason);
    case SuboptAlg::PERSISTENT_LOWMEM:
      return md::stack::SuboptPersistent<true>::IsSupported(cfg, pf, subopt_cfg, reason);
    }
    break;
  }
  if (reason) *reason = fmt::format("not available for {} backend", kind);
  return false;
}

bool PfnAlgIsSupported(BackendKind kind, PfnAlg alg, const erg::EnergyCfg& cfg,
    const erg::PseudofreeCfg& pf, std::string* reason) {
  if (alg == PfnAlg::AUTO) {
    if (reason) *reason = "AUTO is not a concrete algorithm";
    return false;
  }
  if (!BackendIsSupported(kind, cfg, pf, reason)) return false;
  if (alg == PfnAlg::BRUTE) return md::brute::Brute::IsSupported(cfg, pf, reason);

  switch (kind) {
  case BackendKind::BASE:
    switch (alg) {
    case PfnAlg::AUTO:
    case PfnAlg::BRUTE: break;
    case PfnAlg::DEBUG: return md::base::PfnDebug::IsSupported(cfg, pf, reason);
    case PfnAlg::OPT: return md::base::PfnOpt::IsSupported(cfg, pf, reason);
    }
    break;
  case BackendKind::BASEOPT:
    switch (alg) {
    case PfnAlg::AUTO:
    case PfnAlg::BRUTE: break;
    case PfnAlg::DEBUG: return md::base::opt::PfnDebug::IsSupported(cfg, pf, reason);
    case PfnAlg::OPT: return md::base::opt::PfnOpt::IsSupported(cfg, pf, reason);
    }
    break;
  case BackendKind::STACK:
    switch (alg) {
    case PfnAlg::AUTO:
    case PfnAlg::BRUTE:
    case PfnAlg::DEBUG:
    case PfnAlg::OPT: break;
    }
    break;
  }
  if (reason) *reason = fmt::format("not available for {} backend", kind);
  return false;
}

std::optional<MfeAlg> ResolveMfeAlg(
    BackendKind kind, const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf, std::string* log) {
  for (MfeAlg alg : MfePriorityForBackend(kind, /*include_brute=*/false)) {
    std::string reason;
    if (MfeAlgIsSupported(kind, alg, cfg, pf, log ? &reason : nullptr)) {
      return alg;
    }
    if (log) *log += fmt::format("  {}: {}\n", alg, reason);
  }
  return std::nullopt;
}

std::optional<SuboptAlg> ResolveSuboptAlg(BackendKind kind, const erg::EnergyCfg& cfg,
    const erg::PseudofreeCfg& pf, const subopt::SuboptCfg& subopt_cfg, std::string* log) {
  for (SuboptAlg alg : SuboptPriorityForBackend(kind, /*include_brute=*/false)) {
    std::string reason;
    if (SuboptAlgIsSupported(kind, alg, cfg, pf, subopt_cfg, log ? &reason : nullptr)) {
      return alg;
    }
    if (log) *log += fmt::format("  {}: {}\n", alg, reason);
  }
  return std::nullopt;
}

std::optional<PfnAlg> ResolvePfnAlg(
    BackendKind kind, const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf, std::string* log) {
  for (PfnAlg alg : PfnPriorityForBackend(kind, /*include_brute=*/false)) {
    std::string reason;
    if (PfnAlgIsSupported(kind, alg, cfg, pf, log ? &reason : nullptr)) {
      return alg;
    }
    if (log) *log += fmt::format("  {}: {}\n", alg, reason);
  }
  return std::nullopt;
}

void RegisterOptsAlgorithm(ArgParse* args) {
  args->RegisterOpt(OPT_MFE_ALG);
  args->RegisterOpt(OPT_SUBOPT_ALG);
  args->RegisterOpt(OPT_PFN_ALG);
}

}  // namespace mrna
