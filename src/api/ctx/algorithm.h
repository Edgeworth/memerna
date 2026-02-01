// Copyright 2024 Eliot Courtney.
#ifndef API_CTX_ALGORITHM_H_
#define API_CTX_ALGORITHM_H_

#include <optional>
#include <string>

#include "api/ctx/backend_cfg.h"
#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "api/subopt/subopt_cfg.h"
#include "util/argparse.h"
#include "util/container.h"
#include "util/enum.h"

namespace mrna {

MAKE_ENUM(MfeAlg, AUTO, BRUTE, DEBUG, OPT, SPARSE_OPT, LYNGSO_SPARSE_OPT);
MAKE_ENUM(
    SuboptAlg, AUTO, BRUTE, DEBUG, ITERATIVE, ITERATIVE_LOWMEM, PERSISTENT, PERSISTENT_LOWMEM);
MAKE_ENUM(PfnAlg, AUTO, BRUTE, DEBUG, OPT);

[[nodiscard]] bool BackendIsSupported(BackendKind kind, const erg::EnergyCfg& cfg,
    const erg::PseudofreeCfg& pf, std::string* reason = nullptr);

[[nodiscard]] smallvec<MfeAlg, EnumCount<MfeAlg>()> MfePriorityForBackend(BackendKind kind, bool include_brute);
[[nodiscard]] smallvec<SuboptAlg, EnumCount<SuboptAlg>()> SuboptPriorityForBackend(BackendKind kind, bool include_brute);
[[nodiscard]] smallvec<PfnAlg, EnumCount<PfnAlg>()> PfnPriorityForBackend(BackendKind kind, bool include_brute);

[[nodiscard]] bool MfeAlgIsSupported(BackendKind kind, MfeAlg alg, const erg::EnergyCfg& cfg,
    const erg::PseudofreeCfg& pf, std::string* reason = nullptr);
[[nodiscard]] bool SuboptAlgIsSupported(BackendKind kind, SuboptAlg alg, const erg::EnergyCfg& cfg,
    const erg::PseudofreeCfg& pf, const subopt::SuboptCfg& subopt_cfg,
    std::string* reason = nullptr);
[[nodiscard]] bool PfnAlgIsSupported(BackendKind kind, PfnAlg alg, const erg::EnergyCfg& cfg,
    const erg::PseudofreeCfg& pf, std::string* reason = nullptr);

[[nodiscard]] std::optional<MfeAlg> ResolveMfeAlg(BackendKind kind, const erg::EnergyCfg& cfg,
    const erg::PseudofreeCfg& pf, std::string* log = nullptr);
[[nodiscard]] std::optional<SuboptAlg> ResolveSuboptAlg(BackendKind kind, const erg::EnergyCfg& cfg,
    const erg::PseudofreeCfg& pf, const subopt::SuboptCfg& subopt_cfg,
    std::string* log = nullptr);
[[nodiscard]] std::optional<PfnAlg> ResolvePfnAlg(BackendKind kind, const erg::EnergyCfg& cfg,
    const erg::PseudofreeCfg& pf, std::string* log = nullptr);

inline const Opt OPT_MFE_ALG = Opt(Opt::ARG)
                                   .LongName("dp-alg")
                                   .Default(MfeAlg::AUTO)
                                   .ChoiceEnum<MfeAlg>()
                                   .Help("which algorithm for mfe folding");
inline const Opt OPT_SUBOPT_ALG = Opt(Opt::ARG)
                                      .LongName("subopt-alg")
                                      .Default(SuboptAlg::AUTO)
                                      .ChoiceEnum<SuboptAlg>()
                                      .Help("which algorithm for suboptimal folding");
inline const Opt OPT_PFN_ALG = Opt(Opt::ARG)
                                   .LongName("pfn-alg")
                                   .Default(PfnAlg::AUTO)
                                   .ChoiceEnum<PfnAlg>()
                                   .Help("which algorithm for the partition function");

void RegisterOptsAlgorithm(ArgParse* args);

}  // namespace mrna

#endif  // API_CTX_ALGORITHM_H_
