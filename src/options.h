// Copyright 2022 Eliot Courtney.
#ifndef OPTIONS_H_
#define OPTIONS_H_

#include "util/argparse.h"

namespace mrna {

// Some common options
inline const Opt OPT_VERBOSE = Opt().LongName("verbose").ShortName("v").Help("verbose output");
inline const Opt OPT_QUIET = Opt().LongName("quiet").ShortName("q").Help("quiet output");
inline const Opt OPT_AFL = mrna::Opt().LongName("afl").Help("run in afl-fuzz mode");

inline const auto OPT_EFN = mrna::Opt().LongName("efn").ShortName("e").Help("run efn");
inline const auto OPT_MFE = mrna::Opt().LongName("mfe").ShortName("f").Help("run mfe");
inline const auto OPT_SUBOPT = mrna::Opt().LongName("subopt").ShortName("s").Help("run subopt");
inline const auto OPT_PART = mrna::Opt().LongName("part").ShortName("p").Help("run partition");

}  // namespace mrna

#endif  // OPTIONS_H_
