// Copyright 2022 Eliot Courtney.
#ifndef OPTIONS_H_
#define OPTIONS_H_

#include "util/argparse.h"

namespace mrna {

// Some common options
inline const Opt OPT_VERBOSE =
    Opt(Opt::FLAG).LongName("verbose").ShortName("v").Help("verbose output");
inline const Opt OPT_QUIET = Opt(Opt::FLAG).LongName("quiet").ShortName("q").Help("quiet output");
inline const Opt OPT_AFL = Opt(Opt::FLAG).LongName("afl").Help("run in afl-fuzz mode");

inline const Opt OPT_EFN = Opt(Opt::FLAG).LongName("efn").ShortName("e").Help("run efn");
inline const Opt OPT_FOLD = Opt(Opt::FLAG).LongName("fold").ShortName("f").Help("run fold");
inline const Opt OPT_SUBOPT = Opt(Opt::FLAG).LongName("subopt").ShortName("s").Help("run subopt");
inline const Opt OPT_PART = Opt(Opt::FLAG).LongName("part").ShortName("p").Help("run partition");

}  // namespace mrna

#endif  // OPTIONS_H_
