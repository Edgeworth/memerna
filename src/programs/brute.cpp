// Copyright 2022 E.

#include "models/brute/brute.h"

#include <fmt/core.h>

#include <algorithm>
#include <string>

#include "api/brute/brute_cfg.h"
#include "api/energy/model.h"
#include "api/options.h"
#include "api/subopt/subopt_cfg.h"
#include "api/trace/trace.h"
#include "model/ctd.h"
#include "model/secondary.h"
#include "programs/print.h"
#include "util/argparse.h"
#include "util/error.h"

int main(int argc, char* argv[]) {
  mrna::InitProgram();
  mrna::ArgParse args;
  mrna::brute::RegisterOpts(&args);
  args.ParseOrExit(argc, argv);
  const auto em = mrna::erg::FromArgParse(args);
  auto cfg = mrna::brute::BruteCfg::FromArgParse(args);
  if (cfg.mfe) {
    cfg.subopt = true;
    cfg.subopt_cfg.strucs = std::max(cfg.subopt_cfg.strucs, 1);
    cfg.subopt_cfg.sorted = true;
  }

  verify(args.PosSize() == 1, "requires primary sequence");
  auto r = mrna::Primary::FromSeq(args.Pos(0));
  auto res = mrna::md::brute::Brute(r, em, cfg).Run();

  if (args.GetOr(mrna::OPT_FOLD)) {
    const auto& mfe = *res.subopts.begin();
    fmt::print("{}\n", mfe.energy);
    fmt::print("{}\n", mfe.tb.s.ToDb());
    fmt::print("{}\n", mfe.tb.ctd.ToString(mfe.tb.s));
  }

  if (args.GetOr(mrna::OPT_SUBOPT)) {
    for (const auto& s : res.subopts) {
      fmt::print("{}\n", s.energy);
      fmt::print("{}\n", s.tb.s.ToDb());
      fmt::print("{}\n", s.tb.ctd.ToString(s.tb.s));
    }
  }

  if (args.GetOr(mrna::OPT_PART)) {
    fmt::print("q: {}\n", res.part.q);
    fmt::print("p:\n");
    PrintPartition(res.part.p);
    fmt::print("\nprobabilities:\n");
    PrintBoltzProbs(res.part.prob);
  }
}
