// Copyright 2021 Eliot Courtney.
#ifndef API_CTX_CTX_CFG_H_
#define API_CTX_CTX_CFG_H_

#include <iosfwd>
#include <string>
#include <vector>

#include "api/energy/model.h"
#include "util/argparse.h"
#include "util/container.h"

namespace mrna {

inline const Opt OPT_MFE_ALG =
    Opt(Opt::ARG)
        .LongName("dp-alg")
        .Default("auto")
        .Choice({"auto", "opt", "sparse-opt", "lyngso-sparse-opt", "debug", "brute"})
        .Help("which algorithm for mfe folding");
inline const Opt OPT_SUBOPT_ALG = Opt(Opt::ARG)
                                      .LongName("subopt-alg")
                                      .Default("auto")
                                      .Choice({"auto", "iterative", "persistent", "debug", "brute"})
                                      .Help("which algorithm for suboptimal folding");
inline const Opt OPT_PART_ALG = Opt(Opt::ARG)
                                    .LongName("part-alg")
                                    .Default("auto")
                                    .Choice({"auto", "opt", "debug", "brute"})
                                    .Help("which algorithm for the partition function");

void RegisterOpts(ArgParse* args);

struct CtxCfg {
  enum class MfeAlg { DEBUG, OPT, SPARSE_OPT, LYNGSO_SPARSE_OPT, AUTO, BRUTE };
  enum class SuboptAlg { DEBUG, ITERATIVE, PERSISTENT, AUTO, BRUTE };
  enum class PartAlg { DEBUG, OPT, AUTO, BRUTE };

  MfeAlg mfe_alg = MfeAlg::AUTO;
  SuboptAlg subopt_alg = SuboptAlg::AUTO;
  PartAlg part_alg = PartAlg::AUTO;

  [[nodiscard]] static constexpr std::vector<MfeAlg> MfeAlgsForModel(
      const erg::EnergyModelPtr& em) {
    return MfeAlgsForModelKind(erg::Kind(em));
  }

  [[nodiscard]] static constexpr std::vector<SuboptAlg> SuboptAlgsForModel(
      const erg::EnergyModelPtr& em) {
    return SuboptAlgsForModelKind(erg::Kind(em));
  }

  [[nodiscard]] static constexpr std::vector<PartAlg> PartitionAlgsForModel(
      const erg::EnergyModelPtr& em) {
    return PartitionAlgsForModelKind(erg::Kind(em));
  }

  [[nodiscard]] static constexpr std::vector<MfeAlg> MfeAlgsForModelKind(erg::ModelKind kind) {
    switch (kind) {
    case erg::ModelKind::T04_LIKE:
      return {MfeAlg::DEBUG, MfeAlg::OPT, MfeAlg::SPARSE_OPT, MfeAlg::LYNGSO_SPARSE_OPT,
          MfeAlg::AUTO, MfeAlg::BRUTE};
    case erg::ModelKind::T22_LIKE: return {MfeAlg::DEBUG, MfeAlg::BRUTE, MfeAlg::AUTO};
    }
    unreachable();
  }

  [[nodiscard]] static constexpr std::vector<SuboptAlg> SuboptAlgsForModelKind(
      erg::ModelKind kind) {
    switch (kind) {
    case erg::ModelKind::T04_LIKE:
      return {SuboptAlg::DEBUG, SuboptAlg::ITERATIVE, SuboptAlg::PERSISTENT, SuboptAlg::AUTO,
          SuboptAlg::BRUTE};
    case erg::ModelKind::T22_LIKE: return {SuboptAlg::ITERATIVE, SuboptAlg::AUTO, SuboptAlg::BRUTE};
    }
    unreachable();
  }

  [[nodiscard]] static constexpr std::vector<PartAlg> PartitionAlgsForModelKind(
      erg::ModelKind kind) {
    switch (kind) {
    case erg::ModelKind::T04_LIKE:
      return {PartAlg::DEBUG, PartAlg::OPT, PartAlg::AUTO, PartAlg::BRUTE};
    case erg::ModelKind::T22_LIKE: return {PartAlg::AUTO, PartAlg::BRUTE};
    }
    unreachable();
  }

  static CtxCfg FromArgParse(const ArgParse& args);
};

std::istream& operator>>(std::istream& str, CtxCfg::MfeAlg& o);
std::ostream& operator<<(std::ostream& str, const CtxCfg::MfeAlg& o);

std::istream& operator>>(std::istream& str, CtxCfg::SuboptAlg& o);
std::ostream& operator<<(std::ostream& str, const CtxCfg::SuboptAlg& o);

std::istream& operator>>(std::istream& str, CtxCfg::PartAlg& o);
std::ostream& operator<<(std::ostream& str, const CtxCfg::PartAlg& o);

}  // namespace mrna

template <>
struct fmt::formatter<mrna::CtxCfg::MfeAlg> : ostream_formatter {};

template <>
struct fmt::formatter<mrna::CtxCfg::SuboptAlg> : ostream_formatter {};

template <>
struct fmt::formatter<mrna::CtxCfg::PartAlg> : ostream_formatter {};

#endif  // API_CTX_CTX_CFG_H_
