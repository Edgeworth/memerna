// Copyright 2016 Eliot Courtney.
#include <iostream>
#include <memory>

#include "bridge/bridge.h"
#include "energy/load_model.h"
#include "parsing.h"

using mrna::computed_t;
using mrna::energy_t;
using mrna::opt_t;

int main(int argc, char* argv[]) {
  mrna::ArgParse argparse({
      {"v", {"be verbose (if possible)"}},
      {"e", {"run efn"}},
      {"f", {"run fold"}},
      {"subopt-delta", opt_t("maximum energy delta from minimum").Arg("-1")},
  });
  argparse.AddOptions(mrna::bridge::BRIDGE_OPTIONS);
  argparse.AddOptions(mrna::CONTEXT_OPTIONS);
  argparse.AddOptions(mrna::energy::ENERGY_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  verify(argparse.HasFlag("e") + argparse.HasFlag("f") == 1, "require exactly one program flag\n%s",
      argparse.Usage().c_str());

  const auto package = mrna::bridge::RnaPackageFromArgParse(argparse);
  const auto& pos = argparse.GetPositional();
  const bool read_stdin = pos.empty();
  std::deque<std::string> rnaqueue(pos.begin(), pos.end());
  if (argparse.HasFlag("e")) {
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
      const auto secondary = mrna::parsing::ParseDotBracketSecondary(seq, db);
      std::string desc;
      const auto res = package->Efn(secondary, argparse.HasFlag("v") ? &desc : nullptr);
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
      const auto r = mrna::parsing::StringToPrimary(seq);
      int subopt_delta = atoi(argparse.GetOption("subopt-delta").c_str());
      if (subopt_delta >= 0) {
        int num_structures = package->Suboptimal(
            [](const computed_t& c) {
              printf("%d %s\n", c.energy, mrna::parsing::PairsToDotBracket(c.s.p).c_str());
            },
            r, subopt_delta);
        printf("%d suboptimal structures:\n", num_structures);
      } else {
        const auto res = package->Fold(r);
        printf("%d\n%s\n", res.energy, mrna::parsing::PairsToDotBracket(res.s.p).c_str());
      }
    }
  }
}
