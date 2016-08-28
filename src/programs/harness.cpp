#include <iostream>
#include <memory>
#include "parsing.h"
#include "argparse.h"
#include "bridge/bridge.h"
#include "energy/energy_model.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  ArgParse argparse({
      {"v", {"be verbose (if possible)"}},
      {"e", {"run efn"}},
      {"f", {"run fold"}},
  });
  argparse.AddOptions(bridge::BRIDGE_OPTIONS);
  argparse.AddOptions(fold::FOLD_OPTIONS);
  argparse.AddOptions(energy::ENERGY_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  verify_expr(
      argparse.HasFlag("e") + argparse.HasFlag("f") == 1,
      "require exactly one program flag\n%s", argparse.Usage().c_str());
  energy::LoadEnergyModelFromArgParse(argparse);

  auto package = bridge::RnaPackageFromArgParse(argparse);

  const auto& pos = argparse.GetPositional();
  bool read_stdin = pos.empty();
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
      auto frna = parsing::ParseDotBracketRna(seq, db);
      std::string desc;
      auto res = package->Efn(frna, argparse.HasFlag("v") ? &desc : nullptr);
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
      auto rna = parsing::StringToRna(seq);
      auto res = package->Fold(rna);
      printf("%d\n%s\n", res.energy, parsing::PairsToDotBracket(res.p).c_str());
    }
  }
}
