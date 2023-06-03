// Copyright 2016 Eliot Courtney.
#include <fmt/core.h>

#include <deque>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "api/bridge/bridge.h"
#include "api/ctx/ctx.h"
#include "api/energy/energy.h"
#include "api/energy/energy_cfg.h"
#include "api/mfe.h"
#include "api/options.h"
#include "api/part.h"
#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "api/trace/trace.h"
#include "model/energy.h"
#include "model/part.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "programs/print.h"
#include "util/argparse.h"
#include "util/error.h"

int main(int argc, char* argv[]) {
  mrna::InitProgram();
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

  const bool efn = args.GetOr(mrna::OPT_EFN);
  const bool fold = args.GetOr(mrna::OPT_FOLD);
  const bool subopt = args.GetOr(mrna::OPT_SUBOPT);
  const bool part = args.GetOr(mrna::OPT_PART);

  verify(efn + fold + subopt + part == 1, "require exactly one program flag\n{}", args.Usage());
  verify(args.Has(mrna::erg::OPT_SEED) + args.Has(mrna::erg::OPT_MEMERNA_DATA) == 1,
      "require exactly one seed or memerna-data flag\n{}", args.Usage());

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
      fmt::print("{}\n{}", res.energy, desc);
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
              fmt::print("{} {}\n", c.energy, c.tb.s.ToDb());
            },
            r, delta);
        fmt::print("{} suboptimal structures\n", strucs);
      } else if (fold) {
        const auto res = package->Fold(r);
        fmt::print("{}\n{}\n", res.mfe.energy, res.tb.s.ToDb());
      } else if (part) {
        auto res = package->Partition(r);
        fmt::print("q: {}\np:\n", res.part.q);
        PrintPartition(res.part.p);
        fmt::print("\nprobabilities:\n");
        PrintBoltzProbs(res.part.prob);
      }
    }
  }
}
