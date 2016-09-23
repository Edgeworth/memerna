#include <iostream>
#include <memory>
#include "bridge/bridge.h"
#include "energy/load_model.h"
#include "parsing.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  ArgParse argparse({
      {"v", {"be verbose (if possible)"}}, {"e", {"run efn"}}, {"f", {"run fold"}},
      {"subopt-delta", ArgParse::option_t("maximum energy delta from minimum").Arg("-1")},
  });
  argparse.AddOptions(bridge::BRIDGE_OPTIONS);
  argparse.AddOptions(fold::FOLD_OPTIONS);
  argparse.AddOptions(energy::ENERGY_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  verify_expr(argparse.HasFlag("e") + argparse.HasFlag("f") == 1,
      "require exactly one program flag\n%s", argparse.Usage().c_str());

  const auto package = bridge::RnaPackageFromArgParse(argparse);
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
        verify_expr(rnaqueue.size() >= 2, "need even number of positional args");
        seq = rnaqueue.front();
        rnaqueue.pop_front();
        db = rnaqueue.front();
        rnaqueue.pop_front();
      }
      const auto secondary = parsing::ParseDotBracketSecondary(seq, db);
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
      const auto r = parsing::StringToPrimary(seq);
      int subopt_delta = atoi(argparse.GetOption("subopt-delta").c_str());
      if (subopt_delta >= 0) {
        const auto structures = package->Suboptimal(r, subopt_delta);
        printf("%zu suboptimal structures:\n", structures.size());
        for (const auto& subopt : structures) {
          printf("%d %s\n", subopt.energy, parsing::PairsToDotBracket(subopt.s.p).c_str());
        }
      } else {
        const auto res = package->Fold(r);
        printf("%d\n%s\n", res.energy, parsing::PairsToDotBracket(res.s.p).c_str());
      }
    }
  }
}
