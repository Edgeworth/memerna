// Copyright 2024 Eliot Courtney.
#include "api/ctx/algorithm.h"

#include <fmt/core.h>

#include <iterator>
#include <string>
#include <utility>

#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "api/subopt/subopt_cfg.h"
#include "backends/base/energy/boltz_model.h"
#include "backends/base/energy/model.h"
#include "backends/base/mfe/mfe_debug.h"
#include "backends/base/mfe/mfe_exterior.h"
#include "backends/base/mfe/mfe_lyngso_sparse_opt.h"
#include "backends/base/mfe/mfe_opt.h"
#include "backends/base/mfe/mfe_sparse_opt.h"
#include "backends/base/pfn/pfn_debug.h"
#include "backends/base/pfn/pfn_opt.h"
#include "backends/base/subopt/subopt_debug.h"
#include "backends/base/subopt/subopt_iterative.h"
#include "backends/base/subopt/subopt_persistent.h"
#include "backends/base/trace/trace.h"
#include "backends/baseopt/energy/boltz_model.h"
#include "backends/baseopt/energy/model.h"
#include "backends/baseopt/mfe/mfe_debug.h"
#include "backends/baseopt/mfe/mfe_exterior.h"
#include "backends/baseopt/mfe/mfe_lyngso_sparse_opt.h"
#include "backends/baseopt/mfe/mfe_opt.h"
#include "backends/baseopt/mfe/mfe_sparse_opt.h"
#include "backends/baseopt/pfn/pfn_debug.h"
#include "backends/baseopt/pfn/pfn_opt.h"
#include "backends/baseopt/subopt/subopt_debug.h"
#include "backends/baseopt/subopt/subopt_iterative.h"
#include "backends/baseopt/subopt/subopt_persistent.h"
#include "backends/baseopt/trace/trace.h"
#include "backends/common/base/dp.h"
#include "backends/stack/energy/model.h"
#include "backends/stack/mfe/dp.h"
#include "backends/stack/mfe/mfe_exterior.h"
#include "backends/stack/mfe/mfe_opt.h"
#include "backends/stack/subopt/subopt_iterative.h"
#include "backends/stack/subopt/subopt_persistent.h"
#include "backends/stack/trace/trace.h"
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
  case BackendKind::BASE:
    result.assign(std::begin(MFE_PRIORITY_BASE), std::end(MFE_PRIORITY_BASE));
    break;
  case BackendKind::BASEOPT:
    result.assign(std::begin(MFE_PRIORITY_BASEOPT), std::end(MFE_PRIORITY_BASEOPT));
    break;
  case BackendKind::STACK:
    result.assign(std::begin(MFE_PRIORITY_STACK), std::end(MFE_PRIORITY_STACK));
    break;
  }
  if (include_brute) result.push_back(MfeAlg::BRUTE);
  return result;
}

smallvec<SuboptAlg, EnumCount<SuboptAlg>()> SuboptPriorityForBackend(
    BackendKind kind, bool include_brute) {
  smallvec<SuboptAlg, EnumCount<SuboptAlg>()> result;
  switch (kind) {
  case BackendKind::BASE:
    result.assign(std::begin(SUBOPT_PRIORITY_BASE), std::end(SUBOPT_PRIORITY_BASE));
    break;
  case BackendKind::BASEOPT:
    result.assign(std::begin(SUBOPT_PRIORITY_BASEOPT), std::end(SUBOPT_PRIORITY_BASEOPT));
    break;
  case BackendKind::STACK:
    result.assign(std::begin(SUBOPT_PRIORITY_STACK), std::end(SUBOPT_PRIORITY_STACK));
    break;
  }
  if (include_brute) result.push_back(SuboptAlg::BRUTE);
  return result;
}

smallvec<PfnAlg, EnumCount<PfnAlg>()> PfnPriorityForBackend(BackendKind kind, bool include_brute) {
  smallvec<PfnAlg, EnumCount<PfnAlg>()> result;
  switch (kind) {
  case BackendKind::BASE:
    result.assign(std::begin(PFN_PRIORITY_BASE), std::end(PFN_PRIORITY_BASE));
    break;
  case BackendKind::BASEOPT:
    result.assign(std::begin(PFN_PRIORITY_BASEOPT), std::end(PFN_PRIORITY_BASEOPT));
    break;
  case BackendKind::STACK: break;
  }
  if (include_brute) result.push_back(PfnAlg::BRUTE);
  return result;
}

std::optional<MfeFn> GetMfeFn(BackendKind kind, MfeAlg alg) {
  switch (kind) {
  case BackendKind::BASE:
    switch (alg) {
    case MfeAlg::DEBUG:
      return [](const BackendModelPtr& m, const Primary& r, mfe::DpState& dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf) {
        md::base::MfeDebug::Run(
            r, std::get<md::base::Model::Ptr>(m), std::get<md::base::DpState>(dp), cfg, pf);
      };
    case MfeAlg::OPT:
      return [](const BackendModelPtr& m, const Primary& r, mfe::DpState& dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf) {
        md::base::MfeOpt::Run(
            r, std::get<md::base::Model::Ptr>(m), std::get<md::base::DpState>(dp), cfg, pf);
      };
    case MfeAlg::SPARSE_OPT:
      return [](const BackendModelPtr& m, const Primary& r, mfe::DpState& dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf) {
        md::base::MfeSparseOpt::Run(
            r, std::get<md::base::Model::Ptr>(m), std::get<md::base::DpState>(dp), cfg, pf);
      };
    case MfeAlg::LYNGSO_SPARSE_OPT:
      return [](const BackendModelPtr& m, const Primary& r, mfe::DpState& dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf) {
        md::base::MfeLyngsoSparseOpt::Run(
            r, std::get<md::base::Model::Ptr>(m), std::get<md::base::DpState>(dp), cfg, pf);
      };
    case MfeAlg::AUTO:
    case MfeAlg::BRUTE: break;
    }
    break;
  case BackendKind::BASEOPT:
    switch (alg) {
    case MfeAlg::DEBUG:
      return [](const BackendModelPtr& m, const Primary& r, mfe::DpState& dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf) {
        md::base::opt::MfeDebug::Run(
            r, std::get<md::base::opt::Model::Ptr>(m), std::get<md::base::DpState>(dp), cfg, pf);
      };
    case MfeAlg::OPT:
      return [](const BackendModelPtr& m, const Primary& r, mfe::DpState& dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf) {
        md::base::opt::MfeOpt::Run(
            r, std::get<md::base::opt::Model::Ptr>(m), std::get<md::base::DpState>(dp), cfg, pf);
      };
    case MfeAlg::SPARSE_OPT:
      return [](const BackendModelPtr& m, const Primary& r, mfe::DpState& dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf) {
        md::base::opt::MfeSparseOpt::Run(
            r, std::get<md::base::opt::Model::Ptr>(m), std::get<md::base::DpState>(dp), cfg, pf);
      };
    case MfeAlg::LYNGSO_SPARSE_OPT:
      return [](const BackendModelPtr& m, const Primary& r, mfe::DpState& dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf) {
        md::base::opt::MfeLyngsoSparseOpt::Run(
            r, std::get<md::base::opt::Model::Ptr>(m), std::get<md::base::DpState>(dp), cfg, pf);
      };
    case MfeAlg::AUTO:
    case MfeAlg::BRUTE: break;
    }
    break;
  case BackendKind::STACK:
    switch (alg) {
    case MfeAlg::OPT:
      return [](const BackendModelPtr& m, const Primary& r, mfe::DpState& dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf) {
        md::stack::MfeOpt::Run(
            r, std::get<md::stack::Model::Ptr>(m), std::get<md::stack::DpState>(dp), cfg, pf);
      };
    case MfeAlg::AUTO:
    case MfeAlg::BRUTE:
    case MfeAlg::DEBUG:
    case MfeAlg::SPARSE_OPT:
    case MfeAlg::LYNGSO_SPARSE_OPT: break;
    }
    break;
  }
  return std::nullopt;
}

bool MfeAlgIsSupported(BackendKind kind, MfeAlg alg, const erg::EnergyCfg& cfg,
    const erg::PseudofreeCfg& pf, std::string* reason) {
  if (!BackendIsSupported(kind, cfg, pf, reason)) return false;
  if (alg == MfeAlg::AUTO || alg == MfeAlg::BRUTE) return true;

  switch (kind) {
  case BackendKind::BASE:
    switch (alg) {
    case MfeAlg::DEBUG: return md::base::MfeDebug::IsSupported(cfg, pf, reason);
    case MfeAlg::OPT: return md::base::MfeOpt::IsSupported(cfg, pf, reason);
    case MfeAlg::SPARSE_OPT: return md::base::MfeSparseOpt::IsSupported(cfg, pf, reason);
    case MfeAlg::LYNGSO_SPARSE_OPT:
      return md::base::MfeLyngsoSparseOpt::IsSupported(cfg, pf, reason);
    case MfeAlg::AUTO:
    case MfeAlg::BRUTE: break;
    }
    break;
  case BackendKind::BASEOPT:
    switch (alg) {
    case MfeAlg::DEBUG: return md::base::opt::MfeDebug::IsSupported(cfg, pf, reason);
    case MfeAlg::OPT: return md::base::opt::MfeOpt::IsSupported(cfg, pf, reason);
    case MfeAlg::SPARSE_OPT: return md::base::opt::MfeSparseOpt::IsSupported(cfg, pf, reason);
    case MfeAlg::LYNGSO_SPARSE_OPT:
      return md::base::opt::MfeLyngsoSparseOpt::IsSupported(cfg, pf, reason);
    case MfeAlg::AUTO:
    case MfeAlg::BRUTE: break;
    }
    break;
  case BackendKind::STACK:
    switch (alg) {
    case MfeAlg::OPT: return md::stack::MfeOpt::IsSupported(cfg, pf, reason);
    case MfeAlg::AUTO:
    case MfeAlg::BRUTE:
    case MfeAlg::DEBUG:
    case MfeAlg::SPARSE_OPT:
    case MfeAlg::LYNGSO_SPARSE_OPT: break;
    }
    break;
  }
  if (reason) *reason = fmt::format("not available for {} backend", kind);
  return false;
}

MfeExteriorFn GetMfeExteriorFn(BackendKind kind) {
  switch (kind) {
  case BackendKind::BASE:
    return [](const BackendModelPtr& m, const Primary& r, mfe::DpState& dp, erg::EnergyCfg cfg,
               const erg::PseudofreeCfg& pf) {
      auto& state = std::get<md::base::DpState>(dp);
      return md::base::MfeExterior(r, std::get<md::base::Model::Ptr>(m), state, cfg, pf);
    };
  case BackendKind::BASEOPT:
    return [](const BackendModelPtr& m, const Primary& r, mfe::DpState& dp, erg::EnergyCfg cfg,
               const erg::PseudofreeCfg& pf) {
      auto& state = std::get<md::base::DpState>(dp);
      return md::base::opt::MfeExterior(r, std::get<md::base::opt::Model::Ptr>(m), state, cfg, pf);
    };
  case BackendKind::STACK:
    return [](const BackendModelPtr& m, const Primary& r, mfe::DpState& dp, erg::EnergyCfg cfg,
               const erg::PseudofreeCfg& pf) {
      auto& state = std::get<md::stack::DpState>(dp);
      return md::stack::MfeExterior(r, std::get<md::stack::Model::Ptr>(m), state, cfg, pf);
    };
  }
  unreachable();
}

TraceFn GetTraceFn(BackendKind kind) {
  switch (kind) {
  case BackendKind::BASE:
    return [](const BackendModelPtr& m, const Primary& r, const mfe::DpState& dp,
               erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf, const trace::TraceCfg& trace_cfg) {
      const auto& state = std::get<md::base::DpState>(dp);
      return md::base::Traceback(r, std::get<md::base::Model::Ptr>(m), state, cfg, pf, trace_cfg);
    };
  case BackendKind::BASEOPT:
    return [](const BackendModelPtr& m, const Primary& r, const mfe::DpState& dp,
               erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf, const trace::TraceCfg& trace_cfg) {
      const auto& state = std::get<md::base::DpState>(dp);
      return md::base::opt::Traceback(
          r, std::get<md::base::opt::Model::Ptr>(m), state, cfg, pf, trace_cfg);
    };
  case BackendKind::STACK:
    return [](const BackendModelPtr& m, const Primary& r, const mfe::DpState& dp,
               erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf, const trace::TraceCfg& trace_cfg) {
      const auto& state = std::get<md::stack::DpState>(dp);
      return md::stack::Traceback(r, std::get<md::stack::Model::Ptr>(m), state, cfg, pf, trace_cfg);
    };
  }
  unreachable();
}

std::optional<SuboptFn> GetSuboptFn(BackendKind kind, SuboptAlg alg) {
  switch (kind) {
  case BackendKind::BASE:
    switch (alg) {
    case SuboptAlg::DEBUG:
      return [](const BackendModelPtr& m, Primary r, mfe::DpState dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf, const subopt::SuboptCallback& fn,
                 subopt::SuboptCfg subopt_cfg) {
        return md::base::SuboptDebug(std::move(r), std::get<md::base::Model::Ptr>(m),
            std::get<md::base::DpState>(std::move(dp)), cfg, pf, subopt_cfg)
            .Run(fn);
      };
    case SuboptAlg::ITERATIVE:
      return [](const BackendModelPtr& m, Primary r, mfe::DpState dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf, const subopt::SuboptCallback& fn,
                 subopt::SuboptCfg subopt_cfg) {
        return md::base::SuboptIterative<false>(std::move(r), std::get<md::base::Model::Ptr>(m),
            std::get<md::base::DpState>(std::move(dp)), cfg, pf, subopt_cfg)
            .Run(fn);
      };
    case SuboptAlg::ITERATIVE_LOWMEM:
      return [](const BackendModelPtr& m, Primary r, mfe::DpState dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf, const subopt::SuboptCallback& fn,
                 subopt::SuboptCfg subopt_cfg) {
        return md::base::SuboptIterative<true>(std::move(r), std::get<md::base::Model::Ptr>(m),
            std::get<md::base::DpState>(std::move(dp)), cfg, pf, subopt_cfg)
            .Run(fn);
      };
    case SuboptAlg::PERSISTENT:
      return [](const BackendModelPtr& m, Primary r, mfe::DpState dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf, const subopt::SuboptCallback& fn,
                 subopt::SuboptCfg subopt_cfg) {
        return md::base::SuboptPersistent<false>(std::move(r), std::get<md::base::Model::Ptr>(m),
            std::get<md::base::DpState>(std::move(dp)), cfg, pf, subopt_cfg)
            .Run(fn);
      };
    case SuboptAlg::PERSISTENT_LOWMEM:
      return [](const BackendModelPtr& m, Primary r, mfe::DpState dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf, const subopt::SuboptCallback& fn,
                 subopt::SuboptCfg subopt_cfg) {
        return md::base::SuboptPersistent<true>(std::move(r), std::get<md::base::Model::Ptr>(m),
            std::get<md::base::DpState>(std::move(dp)), cfg, pf, subopt_cfg)
            .Run(fn);
      };
    case SuboptAlg::AUTO:
    case SuboptAlg::BRUTE: break;
    }
    break;
  case BackendKind::BASEOPT:
    switch (alg) {
    case SuboptAlg::DEBUG:
      return [](const BackendModelPtr& m, Primary r, mfe::DpState dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf, const subopt::SuboptCallback& fn,
                 subopt::SuboptCfg subopt_cfg) {
        return md::base::opt::SuboptDebug(std::move(r), std::get<md::base::opt::Model::Ptr>(m),
            std::get<md::base::DpState>(std::move(dp)), cfg, pf, subopt_cfg)
            .Run(fn);
      };
    case SuboptAlg::ITERATIVE:
      return [](const BackendModelPtr& m, Primary r, mfe::DpState dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf, const subopt::SuboptCallback& fn,
                 subopt::SuboptCfg subopt_cfg) {
        return md::base::opt::SuboptIterative<false>(std::move(r),
            std::get<md::base::opt::Model::Ptr>(m), std::get<md::base::DpState>(std::move(dp)), cfg,
            pf, subopt_cfg)
            .Run(fn);
      };
    case SuboptAlg::ITERATIVE_LOWMEM:
      return [](const BackendModelPtr& m, Primary r, mfe::DpState dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf, const subopt::SuboptCallback& fn,
                 subopt::SuboptCfg subopt_cfg) {
        return md::base::opt::SuboptIterative<true>(std::move(r),
            std::get<md::base::opt::Model::Ptr>(m), std::get<md::base::DpState>(std::move(dp)), cfg,
            pf, subopt_cfg)
            .Run(fn);
      };
    case SuboptAlg::PERSISTENT:
      return [](const BackendModelPtr& m, Primary r, mfe::DpState dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf, const subopt::SuboptCallback& fn,
                 subopt::SuboptCfg subopt_cfg) {
        return md::base::opt::SuboptPersistent<false>(std::move(r),
            std::get<md::base::opt::Model::Ptr>(m), std::get<md::base::DpState>(std::move(dp)), cfg,
            pf, subopt_cfg)
            .Run(fn);
      };
    case SuboptAlg::PERSISTENT_LOWMEM:
      return [](const BackendModelPtr& m, Primary r, mfe::DpState dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf, const subopt::SuboptCallback& fn,
                 subopt::SuboptCfg subopt_cfg) {
        return md::base::opt::SuboptPersistent<true>(std::move(r),
            std::get<md::base::opt::Model::Ptr>(m), std::get<md::base::DpState>(std::move(dp)), cfg,
            pf, subopt_cfg)
            .Run(fn);
      };
    case SuboptAlg::AUTO:
    case SuboptAlg::BRUTE: break;
    }
    break;
  case BackendKind::STACK:
    switch (alg) {
    case SuboptAlg::ITERATIVE:
      return [](const BackendModelPtr& m, Primary r, mfe::DpState dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf, const subopt::SuboptCallback& fn,
                 subopt::SuboptCfg subopt_cfg) {
        return md::stack::SuboptIterative<false>(std::move(r), std::get<md::stack::Model::Ptr>(m),
            std::get<md::stack::DpState>(std::move(dp)), cfg, pf, subopt_cfg)
            .Run(fn);
      };
    case SuboptAlg::ITERATIVE_LOWMEM:
      return [](const BackendModelPtr& m, Primary r, mfe::DpState dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf, const subopt::SuboptCallback& fn,
                 subopt::SuboptCfg subopt_cfg) {
        return md::stack::SuboptIterative<true>(std::move(r), std::get<md::stack::Model::Ptr>(m),
            std::get<md::stack::DpState>(std::move(dp)), cfg, pf, subopt_cfg)
            .Run(fn);
      };
    case SuboptAlg::PERSISTENT:
      return [](const BackendModelPtr& m, Primary r, mfe::DpState dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf, const subopt::SuboptCallback& fn,
                 subopt::SuboptCfg subopt_cfg) {
        return md::stack::SuboptPersistent<false>(std::move(r), std::get<md::stack::Model::Ptr>(m),
            std::get<md::stack::DpState>(std::move(dp)), cfg, pf, subopt_cfg)
            .Run(fn);
      };
    case SuboptAlg::PERSISTENT_LOWMEM:
      return [](const BackendModelPtr& m, Primary r, mfe::DpState dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf, const subopt::SuboptCallback& fn,
                 subopt::SuboptCfg subopt_cfg) {
        return md::stack::SuboptPersistent<true>(std::move(r), std::get<md::stack::Model::Ptr>(m),
            std::get<md::stack::DpState>(std::move(dp)), cfg, pf, subopt_cfg)
            .Run(fn);
      };
    case SuboptAlg::AUTO:
    case SuboptAlg::BRUTE:
    case SuboptAlg::DEBUG: break;
    }
    break;
  }
  return std::nullopt;
}

bool SuboptAlgIsSupported(BackendKind kind, SuboptAlg alg, const erg::EnergyCfg& cfg,
    const erg::PseudofreeCfg& pf, const subopt::SuboptCfg& subopt_cfg, std::string* reason) {
  if (!BackendIsSupported(kind, cfg, pf, reason)) return false;
  if (alg == SuboptAlg::AUTO || alg == SuboptAlg::BRUTE) return true;

  switch (kind) {
  case BackendKind::BASE:
    switch (alg) {
    case SuboptAlg::DEBUG: return md::base::SuboptDebug::IsSupported(cfg, pf, subopt_cfg, reason);
    case SuboptAlg::ITERATIVE:
      return md::base::SuboptIterative<false>::IsSupported(cfg, pf, subopt_cfg, reason);
    case SuboptAlg::ITERATIVE_LOWMEM:
      return md::base::SuboptIterative<true>::IsSupported(cfg, pf, subopt_cfg, reason);
    case SuboptAlg::PERSISTENT:
      return md::base::SuboptPersistent<false>::IsSupported(cfg, pf, subopt_cfg, reason);
    case SuboptAlg::PERSISTENT_LOWMEM:
      return md::base::SuboptPersistent<true>::IsSupported(cfg, pf, subopt_cfg, reason);
    case SuboptAlg::AUTO:
    case SuboptAlg::BRUTE: break;
    }
    break;
  case BackendKind::BASEOPT:
    switch (alg) {
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
    case SuboptAlg::AUTO:
    case SuboptAlg::BRUTE: break;
    }
    break;
  case BackendKind::STACK:
    switch (alg) {
    case SuboptAlg::ITERATIVE:
      return md::stack::SuboptIterative<false>::IsSupported(cfg, pf, subopt_cfg, reason);
    case SuboptAlg::ITERATIVE_LOWMEM:
      return md::stack::SuboptIterative<true>::IsSupported(cfg, pf, subopt_cfg, reason);
    case SuboptAlg::PERSISTENT:
      return md::stack::SuboptPersistent<false>::IsSupported(cfg, pf, subopt_cfg, reason);
    case SuboptAlg::PERSISTENT_LOWMEM:
      return md::stack::SuboptPersistent<true>::IsSupported(cfg, pf, subopt_cfg, reason);
    case SuboptAlg::AUTO:
    case SuboptAlg::BRUTE:
    case SuboptAlg::DEBUG: break;
    }
    break;
  }
  if (reason) *reason = fmt::format("not available for {} backend", kind);
  return false;
}

std::optional<PfnFn> GetPfnFn(BackendKind kind, PfnAlg alg) {
  switch (kind) {
  case BackendKind::BASE:
    switch (alg) {
    case PfnAlg::DEBUG:
      return [](const BackendModelPtr& m, const Primary& r, pfn::PfnState& dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf) {
        auto& state = std::get<md::base::PfnState>(dp);
        return md::base::PfnDebug::Run(r, std::get<md::base::Model::Ptr>(m), cfg, state, pf);
      };
    case PfnAlg::OPT:
      return [](const BackendModelPtr& m, const Primary& r, pfn::PfnState& dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf) {
        auto& state = std::get<md::base::PfnState>(dp);
        return md::base::PfnOpt::Run(
            r, md::base::BoltzModel::Create(std::get<md::base::Model::Ptr>(m)), cfg, state, pf);
      };
    case PfnAlg::AUTO:
    case PfnAlg::BRUTE: break;
    }
    break;
  case BackendKind::BASEOPT:
    switch (alg) {
    case PfnAlg::DEBUG:
      return [](const BackendModelPtr& m, const Primary& r, pfn::PfnState& dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf) {
        auto& state = std::get<md::base::PfnState>(dp);
        return md::base::opt::PfnDebug::Run(
            r, std::get<md::base::opt::Model::Ptr>(m), cfg, state, pf);
      };
    case PfnAlg::OPT:
      return [](const BackendModelPtr& m, const Primary& r, pfn::PfnState& dp, erg::EnergyCfg cfg,
                 const erg::PseudofreeCfg& pf) {
        auto& state = std::get<md::base::PfnState>(dp);
        return md::base::opt::PfnOpt::Run(r,
            md::base::opt::BoltzModel::Create(std::get<md::base::opt::Model::Ptr>(m)), cfg, state,
            pf);
      };
    case PfnAlg::AUTO:
    case PfnAlg::BRUTE: break;
    }
    break;
  case BackendKind::STACK: break;
  }
  return std::nullopt;
}

bool PfnAlgIsSupported(BackendKind kind, PfnAlg alg, const erg::EnergyCfg& cfg,
    const erg::PseudofreeCfg& pf, std::string* reason) {
  if (!BackendIsSupported(kind, cfg, pf, reason)) return false;
  if (alg == PfnAlg::AUTO || alg == PfnAlg::BRUTE) return true;

  switch (kind) {
  case BackendKind::BASE:
    switch (alg) {
    case PfnAlg::DEBUG: return md::base::PfnDebug::IsSupported(cfg, pf, reason);
    case PfnAlg::OPT: return md::base::PfnOpt::IsSupported(cfg, pf, reason);
    case PfnAlg::AUTO:
    case PfnAlg::BRUTE: break;
    }
    break;
  case BackendKind::BASEOPT:
    switch (alg) {
    case PfnAlg::DEBUG: return md::base::opt::PfnDebug::IsSupported(cfg, pf, reason);
    case PfnAlg::OPT: return md::base::opt::PfnOpt::IsSupported(cfg, pf, reason);
    case PfnAlg::AUTO:
    case PfnAlg::BRUTE: break;
    }
    break;
  case BackendKind::STACK: break;
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
