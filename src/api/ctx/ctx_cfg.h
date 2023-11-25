// Copyright 2021 Eliot Courtney.
#ifndef API_CTX_CTX_CFG_H_
#define API_CTX_CTX_CFG_H_

#include <iosfwd>
#include <string>

#include "api/energy/model.h"
#include "util/argparse.h"
#include "util/container.h"

namespace mrna {

inline const Opt OPT_MFE_ALG = Opt(Opt::ARG)
                                   .LongName("dp-alg")
                                   .Default("fastest")
                                   .Choice({"slowest", "slow", "fastest", "lyngso", "brute"})
                                   .Help("which algorithm for mfe folding");
inline const Opt OPT_SUBOPT_ALG = Opt(Opt::ARG)
                                      .LongName("subopt-alg")
                                      .Default("fastest")
                                      .Choice({"slowest", "fastest", "persistent", "brute"})
                                      .Help("which algorithm for suboptimal folding");
inline const Opt OPT_PART_ALG = Opt(Opt::ARG)
                                    .LongName("part-alg")
                                    .Default("fastest")
                                    .Choice({"slowest", "fastest", "brute"})
                                    .Help("which algorithm for the partition function");

void RegisterOpts(ArgParse* args);

struct CtxCfg {
  enum class MfeAlg { SLOWEST, SLOW, FASTEST, LYNGSO, BRUTE };
  enum class SuboptAlg { SLOWEST, FASTEST, PERSISTENT, BRUTE };
  enum class PartAlg { SLOWEST, FASTEST, BRUTE };

  MfeAlg mfe_alg = MfeAlg::FASTEST;
  SuboptAlg subopt_alg = SuboptAlg::FASTEST;
  PartAlg part_alg = PartAlg::FASTEST;

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
      return {MfeAlg::SLOWEST, MfeAlg::SLOW, MfeAlg::FASTEST, MfeAlg::LYNGSO, MfeAlg::BRUTE};
    case erg::ModelKind::T22_LIKE: return {MfeAlg::SLOWEST, MfeAlg::BRUTE};
    default: bug();
    }
  }

  [[nodiscard]] static constexpr std::vector<SuboptAlg> SuboptAlgsForModelKind(
      erg::ModelKind kind) {
    switch (kind) {
    case erg::ModelKind::T04_LIKE:
      return {SuboptAlg::SLOWEST, SuboptAlg::FASTEST, SuboptAlg::PERSISTENT, SuboptAlg::BRUTE};
    case erg::ModelKind::T22_LIKE: return {SuboptAlg::SLOWEST, SuboptAlg::BRUTE};
    default: bug();
    }
  }

  [[nodiscard]] static constexpr std::vector<PartAlg> PartitionAlgsForModelKind(
      erg::ModelKind kind) {
    switch (kind) {
    case erg::ModelKind::T04_LIKE: return {PartAlg::SLOWEST, PartAlg::FASTEST, PartAlg::BRUTE};
    case erg::ModelKind::T22_LIKE: return {PartAlg::BRUTE};
    default: bug();
    }
  }

  static CtxCfg FromArgParse(const ArgParse& args);
};

std::istream& operator>>(std::istream& str, CtxCfg::MfeAlg& o);
std::istream& operator>>(std::istream& str, CtxCfg::SuboptAlg& o);
std::istream& operator>>(std::istream& str, CtxCfg::PartAlg& o);

}  // namespace mrna

#endif  // API_CTX_CTX_CFG_H_
