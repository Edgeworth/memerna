// Copyright 2016 Eliot Courtney.
#include <iostream>
#include <memory>

#include "bridge/bridge.h"
#include "compute/energy/load_model.h"
#include "model/config.h"
#include "model/parsing.h"

using mrna::computed_t;
using mrna::energy_t;
using mrna::opt_t;

int main(int argc, char* argv[]) {
  mrna::ArgParse args({
      {"v", {"be verbose (if possible)"}},
      {"e", {"run efn"}},
      {"f", {"run fold"}},
      {"subopt-delta", opt_t("maximum energy delta from minimum").Arg("-1")},
  });
  args.AddOptions(mrna::bridge::BRIDGE_OPTS);
  args.AddOptions(mrna::MODEL_OPTS);
  args.AddOptions(mrna::energy::ENERGY_OPTS);
  args.ParseOrExit(argc, argv);
  verify(args.HasFlag("e") + args.HasFlag("f") == 1, "require exactly one program flag\n%s",
      args.Usage().c_str());
  verify(args.HasFlag("seed") + args.HasFlag("memerna-data") == 1,
      "require exactly one seed or memerna-data flag\n%s", args.Usage().c_str());

  const auto package = mrna::bridge::RnaPackageFromArgParse(args);
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
      const auto secondary = mrna::ParseDotBracketSecondary(seq, db);
      std::string desc;
      const auto res = package->Efn(secondary, args.HasFlag("v") ? &desc : nullptr);
      printf("%d\n%s", res, desc.c_str());
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
      const auto r = mrna::StringToPrimary(seq);
      int subopt_delta = atoi(args.GetOption("subopt-delta").c_str());
      if (subopt_delta >= 0) {
        int num_structures = package->Suboptimal(
            [](const computed_t& c) {
              printf("%d %s\n", c.energy, mrna::PairsToDotBracket(c.s.p).c_str());
            },
            r, subopt_delta);
        printf("%d suboptimal structures:\n", num_structures);
      } else {
        const auto res = package->Fold(r);
        printf("%d\n%s\n", res.energy, mrna::PairsToDotBracket(res.s.p).c_str());
      }
    }
  }
}
