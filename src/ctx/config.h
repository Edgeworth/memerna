// Copyright 2021 E.
#ifndef CTX_CONFIG_H_
#define CTX_CONFIG_H_

#include <iosfwd>
#include <string>

#include "util/argparse.h"

namespace mrna::ctx {

inline const Opt OPT_DP_ALG = Opt(Opt::ARG)
                                  .LongName("dp-alg")
                                  .Default("2")
                                  .Choice({"0", "1", "2", "3", "brute"})
                                  .Help("which algorithm for mfe folding");
inline const Opt OPT_SUBOPT_ALG = Opt(Opt::ARG)
                                      .LongName("subopt-alg")
                                      .Default("1")
                                      .Choice({"0", "1", "brute"})
                                      .Help("which algorithm for suboptimal folding");
inline const Opt OPT_PART_ALG = Opt(Opt::ARG)
                                    .LongName("part-alg")
                                    .Default("1")
                                    .Choice({"0", "1", "brute"})
                                    .Help("which algorithm for the partition function");

void RegisterOpts(ArgParse* args);

struct CtxCfg {
  enum class DpAlg { ZERO, ONE, TWO, THREE, BRUTE };
  enum class SuboptAlg { ZERO, ONE, BRUTE };
  enum class PartAlg { ZERO, ONE, BRUTE };

  inline static constexpr DpAlg DP_ALGS[] = {
      DpAlg::ZERO, DpAlg::ONE, DpAlg::TWO, DpAlg::THREE, DpAlg::BRUTE};
  inline static constexpr SuboptAlg SUBOPT_ALGS[] = {
      SuboptAlg::ZERO, SuboptAlg::ONE, SuboptAlg::BRUTE};
  inline static constexpr PartAlg PART_ALGS[] = {PartAlg::ZERO, PartAlg::ONE, PartAlg::BRUTE};

  DpAlg dp_alg = DpAlg::TWO;
  SuboptAlg subopt_alg = SuboptAlg::ONE;
  PartAlg part_alg = PartAlg::ONE;

  static CtxCfg FromArgParse(const ArgParse& args);
};

std::istream& operator>>(std::istream& str, CtxCfg::DpAlg& o);
std::istream& operator>>(std::istream& str, CtxCfg::SuboptAlg& o);
std::istream& operator>>(std::istream& str, CtxCfg::PartAlg& o);

}  // namespace mrna::ctx

#endif  // CTX_CONFIG_H_
