// Copyright 2024 Eliot Courtney.
#ifndef API_CTX_ALGORITHM_H_
#define API_CTX_ALGORITHM_H_

#include <functional>
#include <optional>
#include <string>

#include "api/ctx/backend.h"
#include "api/ctx/backend_cfg.h"
#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "api/mfe.h"
#include "api/pfn.h"
#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "api/trace/trace.h"
#include "api/trace/trace_cfg.h"
#include "util/argparse.h"
#include "util/container.h"
#include "util/enum.h"

namespace mrna {

MAKE_ENUM(MfeAlg, BRUTE, DEBUG, OPT, SPARSE_OPT, LYNGSO_SPARSE_OPT);
MAKE_ENUM(SuboptAlg, BRUTE, DEBUG, ITERATIVE, ITERATIVE_LOWMEM, PERSISTENT, PERSISTENT_LOWMEM);
MAKE_ENUM(PfnAlg, BRUTE, DEBUG, OPT);

struct BackendEfnPriority {
  BackendKind backend;
  int priority;
};

struct BackendMfePriority {
  BackendKind backend;
  MfeAlg alg;
  int priority;
};

struct BackendSuboptPriority {
  BackendKind backend;
  SuboptAlg alg;
  int priority;
};

struct BackendPfnPriority {
  BackendKind backend;
  PfnAlg alg;
  int priority;
};

using MfeFn = std::function<void(const BackendModelPtr&, const Primary&, mfe::DpState&,
    erg::EnergyCfg, const erg::PseudofreeCfg&)>;

using MfeExteriorFn = std::function<Energy(const BackendModelPtr&, const Primary&, mfe::DpState&,
    erg::EnergyCfg, const erg::PseudofreeCfg&)>;

using TraceFn = std::function<trace::TraceResult(const BackendModelPtr&, const Primary&,
    const mfe::DpState&, erg::EnergyCfg, const erg::PseudofreeCfg&, const trace::TraceCfg&)>;

using SuboptFn = std::function<int(const BackendModelPtr&, Primary, mfe::DpState, erg::EnergyCfg,
    const erg::PseudofreeCfg&, const subopt::SuboptCallback&, subopt::SuboptCfg)>;

using PfnFn = std::function<PfnTables(const BackendModelPtr&, const Primary&, pfn::PfnState&,
    erg::EnergyCfg, const erg::PseudofreeCfg&)>;

[[nodiscard]] MfeFn GetMfeFn(BackendKind kind, MfeAlg alg);
[[nodiscard]] MfeExteriorFn GetMfeExteriorFn(BackendKind kind);
[[nodiscard]] TraceFn GetTraceFn(BackendKind kind);
[[nodiscard]] SuboptFn GetSuboptFn(BackendKind kind, SuboptAlg alg);
[[nodiscard]] PfnFn GetPfnFn(BackendKind kind, PfnAlg alg);

[[nodiscard]] bool MfeAlgIsSupported(BackendKind kind, MfeAlg alg, const BackendCfg& backend_cfg,
    const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf, std::string* reason = nullptr);
[[nodiscard]] bool SuboptAlgIsSupported(BackendKind kind, SuboptAlg alg,
    const BackendCfg& backend_cfg, const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf,
    const subopt::SuboptCfg& subopt_cfg, std::string* reason = nullptr);
[[nodiscard]] bool PfnAlgIsSupported(BackendKind kind, PfnAlg alg, const BackendCfg& backend_cfg,
    const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf, std::string* reason = nullptr);

[[nodiscard]] bool BackendIsSupported(BackendKind kind, const BackendCfg& backend_cfg,
    const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf, std::string* reason = nullptr);

[[nodiscard]] smallvec<BackendEfnPriority, static_cast<size_t>(EnumCount<BackendKind>())>
EfnPriorityForBackend(std::optional<BackendKind> kind);
[[nodiscard]] smallvec<BackendMfePriority,
    static_cast<size_t>(EnumCount<BackendKind>() * EnumCount<MfeAlg>())>
MfePriorityForBackend(
    std::optional<BackendKind> kind, std::optional<MfeAlg> alg, bool include_brute);
[[nodiscard]] smallvec<BackendSuboptPriority,
    static_cast<size_t>(EnumCount<BackendKind>() * EnumCount<SuboptAlg>())>
SuboptPriorityForBackend(
    std::optional<BackendKind> kind, std::optional<SuboptAlg> alg, bool include_brute);
[[nodiscard]] smallvec<BackendPfnPriority,
    static_cast<size_t>(EnumCount<BackendKind>() * EnumCount<PfnAlg>())>
PfnPriorityForBackend(
    std::optional<BackendKind> kind, std::optional<PfnAlg> alg, bool include_brute);

[[nodiscard]] std::optional<BackendEfnPriority> ResolveEfn(std::optional<BackendKind> kind,
    const BackendCfg& backend_cfg, const erg::EnergyCfg& cfg, const erg::PseudofreeCfg& pf,
    std::string* log = nullptr);
[[nodiscard]] std::optional<BackendMfePriority> ResolveMfe(std::optional<BackendKind> kind,
    std::optional<MfeAlg> alg, const BackendCfg& backend_cfg, const erg::EnergyCfg& cfg,
    const erg::PseudofreeCfg& pf, std::string* log = nullptr);
[[nodiscard]] std::optional<BackendSuboptPriority> ResolveSubopt(std::optional<BackendKind> kind,
    std::optional<SuboptAlg> alg, const BackendCfg& backend_cfg, const erg::EnergyCfg& cfg,
    const erg::PseudofreeCfg& pf, const subopt::SuboptCfg& subopt_cfg, std::string* log = nullptr);
[[nodiscard]] std::optional<BackendPfnPriority> ResolvePfn(std::optional<BackendKind> kind,
    std::optional<PfnAlg> alg, const BackendCfg& backend_cfg, const erg::EnergyCfg& cfg,
    const erg::PseudofreeCfg& pf, std::string* log = nullptr);

inline const Opt OPT_MFE_ALG = Opt(Opt::ARG).LongName("dp-alg").ChoiceEnum<MfeAlg>().Help(
    "which algorithm for mfe folding (default: auto-select best)");
inline const Opt OPT_SUBOPT_ALG =
    Opt(Opt::ARG)
        .LongName("subopt-alg")
        .ChoiceEnum<SuboptAlg>()
        .Help("which algorithm for suboptimal folding (default: auto-select best)");
inline const Opt OPT_PFN_ALG = Opt(Opt::ARG).LongName("pfn-alg").ChoiceEnum<PfnAlg>().Help(
    "which algorithm for the partition function (default: auto-select best)");

void RegisterOptsAlgorithm(ArgParse* args);

}  // namespace mrna

#endif  // API_CTX_ALGORITHM_H_
