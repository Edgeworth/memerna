// Copyright 2016 Eliot Courtney.
#include <deque>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "bridge/bridge.h"
#include "compute/energy/energy_cfg.h"
#include "compute/energy/model.h"
#include "compute/mfe/mfe.h"
#include "compute/partition/partition.h"
#include "compute/subopt/subopt.h"
#include "compute/subopt/subopt_cfg.h"
#include "compute/traceback/t04/traceback.h"
#include "ctx/ctx.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "options.h"
#include "programs/print.h"
#include "util/argparse.h"
#include "util/error.h"

int main(int argc, char* argv[]) {
  std::ios_base::sync_with_stdio(false);
  mrna::ArgParse args;
  mrna::bridge::RegisterOpts(&args);
  args.RegisterOpt(mrna::OPT_VERBOSE);
  args.RegisterOpt(mrna::OPT_EFN);
  args.RegisterOpt(mrna::OPT_FOLD);
  args.RegisterOpt(mrna::OPT_SUBOPT);
  // Done manually since we only support a subset of options.
  args.RegisterOpt(mrna::subopt::OPT_SUBOPT_DELTA);
  args.RegisterOpt(mrna::OPT_PART);
  args.ParseOrExit(argc, argv);

  bool efn = args.GetOr(mrna::OPT_EFN);
  bool fold = args.GetOr(mrna::OPT_FOLD);
  bool subopt = args.GetOr(mrna::OPT_SUBOPT);
  bool part = args.GetOr(mrna::OPT_PART);

  verify(efn + fold + subopt + part == 1, "require exactly one program flag\n%s",
      args.Usage().c_str());
  verify(args.Has(mrna::energy::OPT_SEED) + args.Has(mrna::energy::OPT_MEMERNA_DATA) == 1,
      "require exactly one seed or memerna-data flag\n%s", args.Usage().c_str());

  auto package = mrna::bridge::RnaPackage::FromArgParse(args);
  std::deque<std::string> q(args.Pos().begin(), args.Pos().end());
  const bool read_stdin = q.empty();
  if (efn) {
    while (true) {
      std::string seq;
      std::string db;
      if (read_stdin) {
        getline(std::cin, seq);
        getline(std::cin, db);
        if (!std::cin) break;
      } else {
        if (q.empty()) break;
        verify(q.size() >= 2, "need even number of positional args");
        seq = q.front();
        q.pop_front();
        db = q.front();
        q.pop_front();
      }
      auto [r, s] = mrna::ParseSeqDb(seq, db);
      std::string desc;
      const auto res = package->Efn(r, s, args.GetOr(mrna::OPT_VERBOSE) ? &desc : nullptr);
      std::cout << res.energy << '\n' << desc;
    }
  } else {
    while (true) {
      std::string seq;
      if (read_stdin) {
        getline(std::cin, seq);
        if (!std::cin) break;
      } else {
        if (q.empty()) break;
        seq = q.front();
        q.pop_front();
      }
      auto r = mrna::Primary::FromSeq(seq);

      if (subopt) {
        auto delta = args.Get<mrna::Energy>(mrna::subopt::OPT_SUBOPT_DELTA);
        int strucs = package->Suboptimal(
            [](const mrna::subopt::SuboptResult& c) {
              std::cout << c.energy << ' ' << c.s.ToDb() << '\n';
            },
            r, delta);
        std::cout << strucs << " suboptimal structures\n";
      } else if (fold) {
        const auto res = package->Fold(r);
        std::cout << res.mfe.energy << '\n' << res.tb.s.ToDb() << '\n';
      } else if (part) {
        auto res = package->Partition(r);
        std::cout << "q: " << res.part.q << '\n';
        std::cout << "p:\n";
        PrintPartition(res.part);
        std::cout << "\nprobabilities:\n";
        PrintBoltzProbs(res.prob);
      }
    }
  }
}
