// Copyright 2024 Eliot Courtney.
#ifndef API_CTX_ALGORITHM_H_
#define API_CTX_ALGORITHM_H_

#include <vector>

#include "api/ctx/backend.h"
#include "api/ctx/backend_cfg.h"
#include "util/argparse.h"
#include "util/enum.h"

namespace mrna {

MAKE_ENUM(MfeAlg, AUTO, BRUTE, DEBUG, OPT, SPARSE_OPT, LYNGSO_SPARSE_OPT);
MAKE_ENUM(
    SuboptAlg, AUTO, BRUTE, DEBUG, ITERATIVE, ITERATIVE_LOWMEM, PERSISTENT, PERSISTENT_LOWMEM);
MAKE_ENUM(PfnAlg, AUTO, BRUTE, DEBUG, OPT);

[[nodiscard]] constexpr std::vector<MfeAlg> MfeAlgsForBackendKind(BackendKind kind) {
  switch (kind) {
  case BackendKind::BASE:
    return {
        MfeAlg::BRUTE,
        MfeAlg::DEBUG,
        MfeAlg::OPT,
        MfeAlg::SPARSE_OPT,
        MfeAlg::LYNGSO_SPARSE_OPT,
    };
  case BackendKind::BASEOPT:
    return {
        MfeAlg::BRUTE,
        MfeAlg::DEBUG,
        MfeAlg::OPT,
        MfeAlg::SPARSE_OPT,
        MfeAlg::LYNGSO_SPARSE_OPT,
    };
  case BackendKind::STACK:
    return {
        MfeAlg::BRUTE,
        MfeAlg::DEBUG,
    };
  }
  unreachable();
}

[[nodiscard]] constexpr std::vector<SuboptAlg> SuboptAlgsForBackendKind(BackendKind kind) {
  switch (kind) {
  case BackendKind::BASE:
    return {
        SuboptAlg::BRUTE,
        SuboptAlg::DEBUG,
        SuboptAlg::ITERATIVE,
        SuboptAlg::ITERATIVE_LOWMEM,
        SuboptAlg::PERSISTENT,
        SuboptAlg::PERSISTENT_LOWMEM,
    };
  case BackendKind::BASEOPT:
    return {
        SuboptAlg::BRUTE,
        SuboptAlg::DEBUG,
        SuboptAlg::ITERATIVE,
        SuboptAlg::ITERATIVE_LOWMEM,
        SuboptAlg::PERSISTENT,
        SuboptAlg::PERSISTENT_LOWMEM,
    };
  case BackendKind::STACK:
    return {
        SuboptAlg::BRUTE,
        SuboptAlg::ITERATIVE,
        SuboptAlg::ITERATIVE_LOWMEM,
        SuboptAlg::PERSISTENT,
        SuboptAlg::PERSISTENT_LOWMEM,
    };
  }
  unreachable();
}

[[nodiscard]] constexpr std::vector<PfnAlg> PfnAlgsForBackendKind(BackendKind kind) {
  switch (kind) {
  case BackendKind::BASE:
    return {
        PfnAlg::BRUTE,
        PfnAlg::DEBUG,
        PfnAlg::OPT,
    };
  case BackendKind::BASEOPT:
    return {
        PfnAlg::BRUTE,
        PfnAlg::DEBUG,
        PfnAlg::OPT,
    };
  case BackendKind::STACK: return {PfnAlg::BRUTE};
  }
  unreachable();
}

// Convenience functions that take BackendModelPtr
[[nodiscard]] inline std::vector<MfeAlg> MfeAlgsForBackend(const BackendModelPtr& m) {
  return MfeAlgsForBackendKind(GetBackendKind(m));
}

[[nodiscard]] inline std::vector<SuboptAlg> SuboptAlgsForBackend(const BackendModelPtr& m) {
  return SuboptAlgsForBackendKind(GetBackendKind(m));
}

[[nodiscard]] inline std::vector<PfnAlg> PfnAlgsForBackend(const BackendModelPtr& m) {
  return PfnAlgsForBackendKind(GetBackendKind(m));
}

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
