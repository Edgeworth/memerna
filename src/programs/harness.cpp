// Copyright 2016 Eliot Courtney.
#include <cstdio>
#include <deque>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "bridge/bridge.h"
#include "compute/energy/config.h"
#include "compute/energy/energy.h"
#include "compute/mfe/mfe.h"
#include "compute/subopt/config.h"
#include "compute/subopt/subopt.h"
#include "compute/traceback/traceback.h"
#include "ctx/ctx.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/argparse.h"
#include "util/error.h"

inline const auto OPT_EFN = mrna::Opt().ShortName("e").Help("run efn");
inline const auto OPT_FOLD = mrna::Opt().ShortName("f").Help("run fold");

int main(int argc, char* argv[]) {
  mrna::ArgParse args;
  mrna::bridge::RegisterOpts(&args);
  args.RegisterOpt(mrna::OPT_VERBOSE);
  args.RegisterOpt(OPT_EFN);
  args.RegisterOpt(OPT_FOLD);
  args.ParseOrExit(argc, argv);

  verify(args.Has(OPT_EFN) + args.Has(OPT_FOLD) == 1, "require exactly one program flag\n%s",
      args.Usage().c_str());
  verify(args.Has(mrna::energy::OPT_SEED) + args.Has(mrna::energy::OPT_MEMERNA_DATA) == 1,
      "require exactly one seed or memerna-data flag\n%s", args.Usage().c_str());

  const auto package = mrna::bridge::RnaPackage::FromArgParse(args);
  std::deque<std::string> q(args.Pos().begin(), args.Pos().end());
  const bool read_stdin = q.empty();
  if (args.Has(OPT_EFN)) {
    while (1) {
      std::string seq, db;
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
      auto [r, s] = mrna::ParsePrimaryDotBracket(seq, db);
      std::string desc;
      const auto res =
          package->Efn(std::move(r), std::move(s), args.Has(mrna::OPT_VERBOSE) ? &desc : nullptr);
      printf("%d\n%s", res.energy, desc.c_str());
    }
  } else {
    while (1) {
      std::string seq;
      if (read_stdin) {
        getline(std::cin, seq);
        if (!std::cin) break;
      } else {
        if (q.empty()) break;
        seq = q.front();
        q.pop_front();
      }
      auto r = mrna::Primary::FromString(seq);
      int subopt_delta = args.Get<int>(mrna::subopt::OPT_SUBOPT_DELTA);
      if (subopt_delta >= 0) {
        int strucs = package->Suboptimal(
            [](const mrna::subopt::SuboptResult& c) {
              printf("%d %s\n", c.energy, c.tb.s.ToDotBracket().c_str());
            },
            std::move(r), subopt_delta);
        printf("%d suboptimal structures:\n", strucs);
      } else {
        const auto res = package->Fold(std::move(r));
        printf("%d\n%s\n", res.mfe.energy, res.tb.s.ToDotBracket().c_str());
      }
    }
  }
}
