// Copyright 2024 Eliot Courtney.
#include "api/ctx/algorithm.h"

#include <fmt/core.h>

#include <algorithm>
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

constexpr int MfePriority(BackendKind backend, MfeAlg alg) {
  switch (backend) {
  case BackendKind::AUTO: return 0;
  case BackendKind::BASE:
    switch (alg) {
    case MfeAlg::AUTO: return 0;
    case MfeAlg::BRUTE: return 912;
    case MfeAlg::DEBUG: return 321;
    case MfeAlg::OPT: return 221;
    case MfeAlg::SPARSE_OPT: return 121;
    case MfeAlg::LYNGSO_SPARSE_OPT: return 122;
    }
    break;
  case BackendKind::BASEOPT:
    switch (alg) {
    case MfeAlg::AUTO: return 0;
    case MfeAlg::BRUTE: return 911;
    case MfeAlg::DEBUG: return 311;
    case MfeAlg::OPT: return 211;
    case MfeAlg::SPARSE_OPT: return 111;
    case MfeAlg::LYNGSO_SPARSE_OPT: return 112;
    }
    break;
  case BackendKind::STACK:
    switch (alg) {
    case MfeAlg::AUTO: return 0;
    case MfeAlg::BRUTE: return 913;
    case MfeAlg::DEBUG: return 0;
    case MfeAlg::OPT: return 231;
    case MfeAlg::SPARSE_OPT: return 0;
    case MfeAlg::LYNGSO_SPARSE_OPT: return 0;
    }
    break;
  }
  unreachable();
}

constexpr int SuboptPriority(BackendKind backend, SuboptAlg alg) {
  switch (backend) {
  case BackendKind::AUTO: return 0;
  case BackendKind::BASE:
    switch (alg) {
    case SuboptAlg::AUTO: return 0;
    case SuboptAlg::BRUTE: return 912;
    case SuboptAlg::DEBUG: return 421;
    case SuboptAlg::ITERATIVE: return 211;
    case SuboptAlg::ITERATIVE_LOWMEM: return 221;
    case SuboptAlg::PERSISTENT: return 212;
    case SuboptAlg::PERSISTENT_LOWMEM: return 222;
    }
    break;
  case BackendKind::BASEOPT:
    switch (alg) {
    case SuboptAlg::AUTO: return 0;
    case SuboptAlg::BRUTE: return 911;
    case SuboptAlg::DEBUG: return 411;
    case SuboptAlg::ITERATIVE: return 111;
    case SuboptAlg::ITERATIVE_LOWMEM: return 121;
    case SuboptAlg::PERSISTENT: return 112;
    case SuboptAlg::PERSISTENT_LOWMEM: return 122;
    }
    break;
  case BackendKind::STACK:
    switch (alg) {
    case SuboptAlg::AUTO: return 0;
    case SuboptAlg::BRUTE: return 913;
    case SuboptAlg::DEBUG: return 0;
    case SuboptAlg::ITERATIVE: return 311;
    case SuboptAlg::ITERATIVE_LOWMEM: return 321;
    case SuboptAlg::PERSISTENT: return 312;
    case SuboptAlg::PERSISTENT_LOWMEM: return 322;
    }
    break;
  }
  unreachable();
}

constexpr int PfnPriority(BackendKind backend, PfnAlg alg) {
  switch (backend) {
  case BackendKind::AUTO: return 0;
  case BackendKind::BASE:
    switch (alg) {
    case PfnAlg::AUTO: return 0;
    case PfnAlg::BRUTE: return 912;
    case PfnAlg::DEBUG: return 221;
    case PfnAlg::OPT: return 121;
    }
    break;
  case BackendKind::BASEOPT:
    switch (alg) {
    case PfnAlg::AUTO: return 0;
    case PfnAlg::BRUTE: return 911;
    case PfnAlg::DEBUG: return 211;
    case PfnAlg::OPT: return 111;
    }
    break;
  case BackendKind::STACK:
    switch (alg) {
    case PfnAlg::AUTO: return 0;
    case PfnAlg::BRUTE: return 913;
    case PfnAlg::DEBUG: return 0;
    case PfnAlg::OPT: return 0;
    }
    break;
  }
  unreachable();
}

}  // namespace

bool BackendIsSupported(BackendKind kind, const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf,
    std::string* reason) {
  switch (kind) {
  case BackendKind::AUTO:
    for (int b = 0; b < EnumCount<BackendKind>(); ++b) {
      auto backend = static_cast<BackendKind>(b);
      if (backend == BackendKind::AUTO) continue;
      if (BackendIsSupported(backend, cfg, pf, nullptr)) return true;
    }
    if (reason) *reason = "no backend supports this configuration";
    return false;
  case BackendKind::BASE: return md::base::Model::IsSupported(cfg, pf, reason);
  case BackendKind::BASEOPT: return md::base::opt::Model::IsSupported(cfg, pf, reason);
  case BackendKind::STACK: return md::stack::Model::IsSupported(cfg, pf, reason);
  }
  unreachable();
}

smallvec<BackendMfePriority, static_cast<size_t>(EnumCount<BackendKind>() * EnumCount<MfeAlg>())>
MfePriorityForBackend(BackendKind kind, bool include_brute) {
  smallvec<BackendMfePriority, static_cast<size_t>(EnumCount<BackendKind>() * EnumCount<MfeAlg>())>
      result;

  for (int b = 0; b < EnumCount<BackendKind>(); ++b) {
    auto backend = static_cast<BackendKind>(b);
    if (backend == BackendKind::AUTO) continue;
    if (kind != BackendKind::AUTO && backend != kind) continue;

    for (int a = 0; a < EnumCount<MfeAlg>(); ++a) {
      auto alg = static_cast<MfeAlg>(a);
      if (alg == MfeAlg::AUTO) continue;
      if (alg == MfeAlg::BRUTE && !include_brute) continue;
      int priority = MfePriority(backend, alg);
      if (priority > 0) result.push_back({.backend = backend, .alg = alg, .priority = priority});
    }
  }
  std::sort(result.begin(), result.end(),
      [](const auto& a, const auto& b) {
        if (a.priority != b.priority) return a.priority < b.priority;
        if (a.backend != b.backend)
          return static_cast<int>(a.backend) < static_cast<int>(b.backend);
        return static_cast<int>(a.alg) < static_cast<int>(b.alg);
      });
  return result;
}

smallvec<BackendSuboptPriority,
    static_cast<size_t>(EnumCount<BackendKind>() * EnumCount<SuboptAlg>())>
SuboptPriorityForBackend(BackendKind kind, bool include_brute) {
  smallvec<BackendSuboptPriority,
      static_cast<size_t>(EnumCount<BackendKind>() * EnumCount<SuboptAlg>())>
      result;

  for (int b = 0; b < EnumCount<BackendKind>(); ++b) {
    auto backend = static_cast<BackendKind>(b);
    if (backend == BackendKind::AUTO) continue;
    if (kind != BackendKind::AUTO && backend != kind) continue;

    for (int a = 0; a < EnumCount<SuboptAlg>(); ++a) {
      auto alg = static_cast<SuboptAlg>(a);
      if (alg == SuboptAlg::AUTO) continue;
      if (alg == SuboptAlg::BRUTE && !include_brute) continue;
      int priority = SuboptPriority(backend, alg);
      if (priority > 0) result.push_back({.backend = backend, .alg = alg, .priority = priority});
    }
  }
  std::sort(result.begin(), result.end(),
      [](const auto& a, const auto& b) {
        if (a.priority != b.priority) return a.priority < b.priority;
        if (a.backend != b.backend)
          return static_cast<int>(a.backend) < static_cast<int>(b.backend);
        return static_cast<int>(a.alg) < static_cast<int>(b.alg);
      });
  return result;
}

smallvec<BackendPfnPriority, static_cast<size_t>(EnumCount<BackendKind>() * EnumCount<PfnAlg>())>
PfnPriorityForBackend(BackendKind kind, bool include_brute) {
  smallvec<BackendPfnPriority, static_cast<size_t>(EnumCount<BackendKind>() * EnumCount<PfnAlg>())>
      result;

  for (int b = 0; b < EnumCount<BackendKind>(); ++b) {
    auto backend = static_cast<BackendKind>(b);
    if (backend == BackendKind::AUTO) continue;
    if (kind != BackendKind::AUTO && backend != kind) continue;

    for (int a = 0; a < EnumCount<PfnAlg>(); ++a) {
      auto alg = static_cast<PfnAlg>(a);
      if (alg == PfnAlg::AUTO) continue;
      if (alg == PfnAlg::BRUTE && !include_brute) continue;
      int priority = PfnPriority(backend, alg);
      if (priority > 0) result.push_back({.backend = backend, .alg = alg, .priority = priority});
    }
  }
  std::sort(result.begin(), result.end(),
      [](const auto& a, const auto& b) {
        if (a.priority != b.priority) return a.priority < b.priority;
        if (a.backend != b.backend)
          return static_cast<int>(a.backend) < static_cast<int>(b.backend);
        return static_cast<int>(a.alg) < static_cast<int>(b.alg);
      });
  return result;
}

std::optional<MfeFn> GetMfeFn(BackendKind kind, MfeAlg alg) {
  verify(kind != BackendKind::AUTO, "must resolve AUTO backend before calling GetMfeFn");
  switch (kind) {
  case BackendKind::AUTO: unreachable();
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
  verify(kind != BackendKind::AUTO, "must resolve AUTO backend before calling MfeAlgIsSupported");
  if (!BackendIsSupported(kind, cfg, pf, reason)) return false;
  if (alg == MfeAlg::BRUTE) return true;
  if (alg == MfeAlg::AUTO) return ResolveMfeAlg(kind, cfg, pf, nullptr).has_value();

  switch (kind) {
  case BackendKind::AUTO: unreachable();
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
  verify(kind != BackendKind::AUTO, "must resolve AUTO backend before calling GetMfeExteriorFn");
  switch (kind) {
  case BackendKind::AUTO: unreachable();
  case BackendKind::BASE:
    return [](const BackendModelPtr& m, const Primary& r, mfe::DpState& dp, erg::EnergyCfg cfg,
               const erg::PseudofreeCfg& pf) {
      auto& state = std::get<md::base::DpState>(dp);
      return md::base::MfeExterior(r, std::get<md::base::Model::Ptr>(m), state, cfg, pf);
    };
  case BackendKind::BASEOPT:
    return [](const BackendModelPtr& m, const Primary& r, mfe::DpState& dp, erg::EnergyCfg /*cfg*/,
               const erg::PseudofreeCfg& pf) {
      auto& state = std::get<md::base::DpState>(dp);
      return md::base::opt::MfeExterior(r, std::get<md::base::opt::Model::Ptr>(m), state, pf);
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
  verify(kind != BackendKind::AUTO, "must resolve AUTO backend before calling GetTraceFn");
  switch (kind) {
  case BackendKind::AUTO: unreachable();
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
  verify(kind != BackendKind::AUTO, "must resolve AUTO backend before calling GetSuboptFn");
  switch (kind) {
  case BackendKind::AUTO: unreachable();
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
  verify(
      kind != BackendKind::AUTO, "must resolve AUTO backend before calling SuboptAlgIsSupported");
  if (!BackendIsSupported(kind, cfg, pf, reason)) return false;
  if (alg == SuboptAlg::BRUTE) return true;
  if (alg == SuboptAlg::AUTO)
    return ResolveSuboptAlg(kind, cfg, pf, subopt_cfg, nullptr).has_value();

  switch (kind) {
  case BackendKind::AUTO: unreachable();
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
  verify(kind != BackendKind::AUTO, "must resolve AUTO backend before calling GetPfnFn");
  switch (kind) {
  case BackendKind::AUTO: unreachable();
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
  verify(kind != BackendKind::AUTO, "must resolve AUTO backend before calling PfnAlgIsSupported");
  if (!BackendIsSupported(kind, cfg, pf, reason)) return false;
  if (alg == PfnAlg::BRUTE) return true;
  if (alg == PfnAlg::AUTO) return ResolvePfnAlg(kind, cfg, pf, nullptr).has_value();

  switch (kind) {
  case BackendKind::AUTO: unreachable();
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

std::optional<BackendMfePriority> ResolveMfeAlg(
    BackendKind kind, const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf, std::string* log) {
  for (const auto& entry : MfePriorityForBackend(kind, /*include_brute=*/false)) {
    std::string reason;
    if (!MfeAlgIsSupported(entry.backend, entry.alg, cfg, pf, log ? &reason : nullptr)) {
      if (log) *log += fmt::format("  {}/{}: {}\n", entry.backend, entry.alg, reason);
      continue;
    }
    return entry;
  }
  return std::nullopt;
}

std::optional<BackendSuboptPriority> ResolveSuboptAlg(BackendKind kind, const erg::EnergyCfg& cfg,
    const erg::PseudofreeCfg& pf, const subopt::SuboptCfg& subopt_cfg, std::string* log) {
  for (const auto& entry : SuboptPriorityForBackend(kind, /*include_brute=*/false)) {
    std::string reason;
    if (!SuboptAlgIsSupported(
            entry.backend, entry.alg, cfg, pf, subopt_cfg, log ? &reason : nullptr)) {
      if (log) *log += fmt::format("  {}/{}: {}\n", entry.backend, entry.alg, reason);
      continue;
    }
    return entry;
  }
  return std::nullopt;
}

std::optional<BackendPfnPriority> ResolvePfnAlg(
    BackendKind kind, const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf, std::string* log) {
  for (const auto& entry : PfnPriorityForBackend(kind, /*include_brute=*/false)) {
    std::string reason;
    if (!PfnAlgIsSupported(entry.backend, entry.alg, cfg, pf, log ? &reason : nullptr)) {
      if (log) *log += fmt::format("  {}/{}: {}\n", entry.backend, entry.alg, reason);
      continue;
    }
    return entry;
  }
  return std::nullopt;
}

void RegisterOptsAlgorithm(ArgParse* args) {
  args->RegisterOpt(OPT_MFE_ALG);
  args->RegisterOpt(OPT_SUBOPT_ALG);
  args->RegisterOpt(OPT_PFN_ALG);
}

}  // namespace mrna
