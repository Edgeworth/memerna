// Copyright 2022 E.
#ifndef MODEL_CONFIG_H_
#define MODEL_CONFIG_H_

#include "util/argparse.h"

namespace mrna {

inline const auto OPT_LONELY_PAIRS =
    mrna::Opt().LongName("lonely-pairs").Help("allow lonely pairs");
inline const auto OPT_CTD = mrna::Opt()
                                .LongName("ctd")
                                .Choice({"none", "no_coax", "all"})
                                .Default("all")
                                .Help("whether to use CTDs");

void RegisterOpts(ArgParse* args);

struct ModelCfg {
  // TODO: Implement and use these options.
  enum class Ctd {
    NONE,  //  Do not use CTDs in folding, subopt, partition, etc.
    NO_COAX,  //  Use only terminal mismatches and dangling ends in folding, subopt, partition, etc.
    ALL,  //  Use CTDs in folding, subopt, partition, etc.
  };

  bool lonely_pairs = false;  // Whether to allow lonely pairs in folding, subopt, partition, etc.
  Ctd ctd = Ctd::ALL;  // Whether to use CTDs in folding, subopt, partition, etc.

  static ModelCfg FromArgParse(const ArgParse& args);
};

std::istream& operator>>(std::istream& str, ModelCfg::Ctd& o);

}  // namespace mrna

#endif  // MODEL_CONFIG_H_
