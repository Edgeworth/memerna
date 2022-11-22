// Copyright 2021 Eliot Courtney.

#ifndef CTX_CTX_CFG_H_
#define CTX_CTX_CFG_H_
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
  enum class DpAlg { SLOWEST, SLOW, FASTEST, LYNGSO, BRUTE };
  enum class SuboptAlg { SLOWEST, FASTEST, BRUTE };
  enum class PartAlg { SLOWEST, FASTEST, BRUTE };

  inline static constexpr DpAlg DP_ALGS[] = {
      DpAlg::SLOWEST, DpAlg::SLOW, DpAlg::FASTEST, DpAlg::LYNGSO, DpAlg::BRUTE};
  inline static constexpr SuboptAlg SUBOPT_ALGS[] = {
      SuboptAlg::SLOWEST, SuboptAlg::FASTEST, SuboptAlg::BRUTE};
  inline static constexpr PartAlg PART_ALGS[] = {
      PartAlg::SLOWEST, PartAlg::FASTEST, PartAlg::BRUTE};

  DpAlg dp_alg = DpAlg::FASTEST;
  SuboptAlg subopt_alg = SuboptAlg::FASTEST;
  PartAlg part_alg = PartAlg::FASTEST;

  static CtxCfg FromArgParse(const ArgParse& args);
};

std::istream& operator>>(std::istream& str, CtxCfg::DpAlg& o);
std::istream& operator>>(std::istream& str, CtxCfg::SuboptAlg& o);
std::istream& operator>>(std::istream& str, CtxCfg::PartAlg& o);

}  // namespace mrna::ctx
#endif  // CTX_CTX_CFG_H_
