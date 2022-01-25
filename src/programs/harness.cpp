// Copyright 2016 E.
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "bridge/bridge.h"
#include "compute/energy/energy.h"
#include "compute/mfe/mfe.h"
#include "compute/subopt/subopt.h"
#include "compute/traceback/traceback.h"
#include "model/config.h"
#include "model/context.h"
#include "model/model.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/argparse.h"
#include "util/error.h"

using mrna::Energy;
using mrna::Opt;

int main(int argc, char* argv[]) {
  mrna::ArgParse args({
      {"v", {"be verbose (if possible)"}},
      {"e", {"run efn"}},
      {"f", {"run fold"}},
      {"subopt-delta", Opt("maximum energy delta from minimum").Arg("-1")},
  });
  args.AddOptions(mrna::bridge::BRIDGE_OPTS);
  args.AddOptions(mrna::MODEL_OPTS);
  args.AddOptions(mrna::energy::ENERGY_OPTS);
  args.ParseOrExit(argc, argv);
  verify(args.HasFlag("e") + args.HasFlag("f") == 1, "require exactly one program flag\n%s",
      args.Usage().c_str());
  verify(args.HasFlag("seed") + args.HasFlag("memerna-data") == 1,
      "require exactly one seed or memerna-data flag\n%s", args.Usage().c_str());

  const auto package = mrna::bridge::RnaPackage::FromArgParse(args);
  const auto& pos = args.GetPositional();
  const bool read_stdin = pos.empty();
  std::deque<std::string> rnaqueue(pos.begin(), pos.end());
  if (args.HasFlag("e")) {
    while (1) {
      std::string seq, db;
      if (read_stdin) {
        getline(std::cin, seq);
        getline(std::cin, db);
        if (!std::cin) break;
      } else {
        if (rnaqueue.empty()) break;
        verify(rnaqueue.size() >= 2, "need even number of positional args");
        seq = rnaqueue.front();
        rnaqueue.pop_front();
        db = rnaqueue.front();
        rnaqueue.pop_front();
      }
      auto [r, s] = mrna::ParsePrimaryDotBracket(seq, db);
      std::string desc;
      const auto res =
          package->Efn(std::move(r), std::move(s), args.HasFlag("v") ? &desc : nullptr);
      printf("%d\n%s", res.energy, desc.c_str());
    }
  } else {
    while (1) {
      std::string seq;
      if (read_stdin) {
        getline(std::cin, seq);
        if (!std::cin) break;
      } else {
        if (rnaqueue.empty()) break;
        seq = rnaqueue.front();
        rnaqueue.pop_front();
      }
      auto r = mrna::Primary::FromString(seq);
      int subopt_delta = atoi(args.GetOption("subopt-delta").c_str());
      if (subopt_delta >= 0) {
        int num_structures = package->Suboptimal(
            [](const mrna::subopt::SuboptResult& c) {
              printf("%d %s\n", c.energy, c.tb.s.ToDotBracket().c_str());
            },
            std::move(r), subopt_delta);
        printf("%d suboptimal structures:\n", num_structures);
      } else {
        const auto res = package->Fold(std::move(r));
        printf("%d\n%s\n", res.mfe.energy, res.tb.s.ToDotBracket().c_str());
      }
    }
  }
}
