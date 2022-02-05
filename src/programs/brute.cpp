// Copyright 2022 Eliot Courtney.
#include "compute/brute/brute.h"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>

#include "compute/brute/config.h"
#include "compute/energy/model.h"
#include "compute/subopt/config.h"
#include "compute/traceback/traceback.h"
#include "model/ctd.h"
#include "model/secondary.h"
#include "options.h"
#include "programs/print.h"
#include "util/argparse.h"
#include "util/error.h"

int main(int argc, char* argv[]) {
  mrna::ArgParse args;
  mrna::brute::RegisterOpts(&args);
  args.ParseOrExit(argc, argv);
  const auto em = mrna::energy::EnergyModel::FromArgParse(args);
  auto cfg = mrna::brute::BruteCfg::FromArgParse(args);
  if (cfg.mfe) {
    cfg.subopt = true;
    cfg.subopt_cfg.strucs = std::max(cfg.subopt_cfg.strucs, 1);
    cfg.subopt_cfg.sorted = true;
  }

  verify(args.PosSize() == 1, "requires primary sequence");
  auto r = mrna::Primary::FromSeq(args.Pos(0));

  auto res = mrna::brute::BruteForce(r, em, std::move(cfg)).Run();

  if (args.GetOr(mrna::OPT_FOLD)) {
    const auto& mfe = *res.subopts.begin();
    printf("%d\n", mfe.energy);
    puts(mfe.tb.s.ToDb().c_str());
    puts(mfe.tb.ctd.ToString(mfe.tb.s).c_str());
  }

  if (args.GetOr(mrna::OPT_SUBOPT)) {
    for (const auto& s : res.subopts) {
      printf("%d\n", s.energy);
      puts(s.tb.s.ToDb().c_str());
      puts(s.tb.ctd.ToString(s.tb.s).c_str());
    }
  }

  if (args.GetOr(mrna::OPT_PART)) {
    std::cout << "q: " << res.part.q << '\n';
    std::cout << "p:\n";
    PrintPartition(res.part);
    std::cout << "\nprobabilities:\n";
    PrintBoltzProbs(res.prob);
  }
}
