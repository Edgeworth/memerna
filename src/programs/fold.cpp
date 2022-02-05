// Copyright 2016 E.
#include <cstdio>
#include <string>

#include "compute/mfe/mfe.h"
#include "compute/traceback/traceback.h"
#include "ctx/config.h"
#include "ctx/ctx.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/argparse.h"
#include "util/error.h"

int main(int argc, char* argv[]) {
  mrna::ArgParse args;
  mrna::ctx::RegisterOpts(&args);
  args.ParseOrExit(argc, argv);
  verify(args.PosSize() == 1, "need primary sequence to fold");

  auto ctx = mrna::ctx::Ctx::FromArgParse(args);
  const auto res = ctx.Fold(mrna::Primary::FromSeq(args.Pos(0)));

  printf("%d\n%s\n%s\n", res.mfe.energy, res.tb.s.ToDb().c_str(),
      res.tb.ctd.ToString(res.tb.s).c_str());
}
