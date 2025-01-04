// Copyright 2021 Eliot Courtney.
#ifndef API_CTX_CTX_CFG_H_
#define API_CTX_CTX_CFG_H_

#include <string>
#include <vector>

#include "api/ctx/backend.h"
#include "util/argparse.h"

namespace mrna {

void RegisterOpts(ArgParse* args);

struct CtxCfg {
  MAKE_NESTED_ENUM(MfeAlg, AUTO, BRUTE, DEBUG, OPT, SPARSE_OPT, LYNGSO_SPARSE_OPT);
  MAKE_NESTED_ENUM(SuboptAlg, AUTO, BRUTE, DEBUG, ITERATIVE, PERSISTENT);
  MAKE_NESTED_ENUM(PfnAlg, AUTO, BRUTE, DEBUG, OPT);

  MfeAlg mfe_alg = MfeAlg::AUTO;

  SuboptAlg subopt_alg = SuboptAlg::AUTO;

  PfnAlg pfn_alg = PfnAlg::AUTO;

  [[nodiscard]] static constexpr std::vector<MfeAlg> MfeAlgsForBackend(const BackendModelPtr& m) {
    return MfeAlgsForBackendKind(GetBackendKind(m));
  }

  [[nodiscard]] static constexpr std::vector<SuboptAlg> SuboptAlgsForBackend(
      const BackendModelPtr& m) {
    return SuboptAlgsForBackendKind(GetBackendKind(m));
  }

  [[nodiscard]] static constexpr std::vector<PfnAlg> PfnAlgsForBackend(const BackendModelPtr& m) {
    return PfnAlgsForBackendKind(GetBackendKind(m));
  }

  [[nodiscard]] static constexpr std::vector<MfeAlg> MfeAlgsForBackendKind(BackendKind kind) {
    switch (kind) {
    case BackendKind::BASE:
      return {
          MfeAlg::AUTO,
          MfeAlg::BRUTE,
          MfeAlg::DEBUG,
          MfeAlg::OPT,
          MfeAlg::SPARSE_OPT,
          MfeAlg::LYNGSO_SPARSE_OPT,
      };
    case BackendKind::BASEOPT:
      return {
          MfeAlg::AUTO,
          MfeAlg::BRUTE,
          MfeAlg::DEBUG,
          MfeAlg::OPT,
          MfeAlg::SPARSE_OPT,
          MfeAlg::LYNGSO_SPARSE_OPT,
      };
    case BackendKind::STACK:
      return {
          MfeAlg::AUTO,
          MfeAlg::BRUTE,
          MfeAlg::DEBUG,
      };
    }
    unreachable();
  }

  [[nodiscard]] static constexpr std::vector<SuboptAlg> SuboptAlgsForBackendKind(BackendKind kind) {
    switch (kind) {
    case BackendKind::BASE:
      return {
          SuboptAlg::AUTO,
          SuboptAlg::BRUTE,
          SuboptAlg::DEBUG,
          SuboptAlg::ITERATIVE,
          SuboptAlg::PERSISTENT,
      };
    case BackendKind::BASEOPT:
      return {
          SuboptAlg::AUTO,
          SuboptAlg::BRUTE,
          SuboptAlg::DEBUG,
          SuboptAlg::ITERATIVE,
          SuboptAlg::PERSISTENT,
      };
    case BackendKind::STACK:
      return {
          SuboptAlg::AUTO,
          SuboptAlg::BRUTE,
          SuboptAlg::ITERATIVE,
          SuboptAlg::PERSISTENT,
      };
    }
    unreachable();
  }

  [[nodiscard]] static constexpr std::vector<PfnAlg> PfnAlgsForBackendKind(BackendKind kind) {
    switch (kind) {
    case BackendKind::BASE:
      return {
          PfnAlg::AUTO,
          PfnAlg::BRUTE,
          PfnAlg::DEBUG,
          PfnAlg::OPT,
      };
    case BackendKind::BASEOPT:
      return {
          PfnAlg::AUTO,
          PfnAlg::BRUTE,
          PfnAlg::DEBUG,
          PfnAlg::OPT,
      };
    case BackendKind::STACK: return {PfnAlg::BRUTE};
    }
    unreachable();
  }

  static CtxCfg FromArgParse(const ArgParse& args);
};

inline const Opt OPT_MFE_ALG = Opt(Opt::ARG)
                                   .LongName("dp-alg")
                                   .Default(CtxCfg::MfeAlg::AUTO)
                                   .ChoiceEnum<CtxCfg::MfeAlg>()
                                   .Help("which algorithm for mfe folding");
inline const Opt OPT_SUBOPT_ALG = Opt(Opt::ARG)
                                      .LongName("subopt-alg")
                                      .Default(CtxCfg::SuboptAlg::AUTO)
                                      .ChoiceEnum<CtxCfg::SuboptAlg>()
                                      .Help("which algorithm for suboptimal folding");
inline const Opt OPT_PFN_ALG = Opt(Opt::ARG)
                                   .LongName("pfn-alg")
                                   .Default(CtxCfg::PfnAlg::AUTO)
                                   .ChoiceEnum<CtxCfg::PfnAlg>()
                                   .Help("which algorithm for the partition function");

}  // namespace mrna

#endif  // API_CTX_CTX_CFG_H_
