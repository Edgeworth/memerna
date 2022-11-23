// Copyright 2022 Eliot Courtney.
#include "compute/brute/brute.h"

#include <algorithm>
#include <iostream>
#include <string>

#include "compute/brute/brute_cfg.h"
#include "compute/energy/model.h"
#include "compute/subopt/subopt_cfg.h"
#include "compute/traceback/t04/traceback.h"
#include "model/ctd.h"
#include "model/secondary.h"
#include "options.h"
#include "programs/print.h"
#include "util/argparse.h"
#include "util/error.h"

int main(int argc, char* argv[]) {
  std::ios_base::sync_with_stdio(false);
  mrna::ArgParse args;
  mrna::brute::RegisterOpts(&args);
  args.ParseOrExit(argc, argv);
  const auto em = mrna::energy::FromArgParse(args);
  auto cfg = mrna::brute::BruteCfg::FromArgParse(args);
  if (cfg.mfe) {
    cfg.subopt = true;
    cfg.subopt_cfg.strucs = std::max(cfg.subopt_cfg.strucs, 1);
    cfg.subopt_cfg.sorted = true;
  }

  verify(args.PosSize() == 1, "requires primary sequence");
  auto r = mrna::Primary::FromSeq(args.Pos(0));
  auto res = mrna::brute::Brute(r, em, cfg).Run();

  if (args.GetOr(mrna::OPT_FOLD)) {
    const auto& mfe = *res.subopts.begin();
    std::cout << mfe.energy << '\n';
    std::cout << mfe.tb.s.ToDb() << '\n';
    std::cout << mfe.tb.ctd.ToString(mfe.tb.s) << '\n';
  }

  if (args.GetOr(mrna::OPT_SUBOPT)) {
    for (const auto& s : res.subopts) {
      std::cout << s.energy << '\n';
      std::cout << s.tb.s.ToDb() << '\n';
      std::cout << s.tb.ctd.ToString(s.tb.s) << '\n';
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
