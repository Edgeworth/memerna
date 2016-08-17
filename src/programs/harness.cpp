#include <iostream>
#include <memory>
#include "parsing.h"
#include "argparse.h"
#include "bridge/bridge.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  ArgParse argparse({
                        {"v", {"be verbose (if possible)"}},
                        {"r", {"rnastructure"}},
                        {"m", {"rnark"}},
                        {"k", {"memerna"}},
                        {"e", {"run efn"}},
                        {"f", {"run fold"}},
                    });

  auto ret = argparse.Parse(argc, argv);

  verify_expr(
      ret.size() == 0,
      "%s\n%s\n", ret.c_str(), argparse.Usage().c_str());
  verify_expr(
      argparse.HasFlag("r") + argparse.HasFlag("m") + argparse.HasFlag("k") == 1,
      "require exactly one package flag\n%s", argparse.Usage().c_str());
  verify_expr(
      argparse.HasFlag("e") + argparse.HasFlag("f") == 1,
      "require exactly one program flag\n%s", argparse.Usage().c_str());

  std::unique_ptr<bridge::RnaPackage> package;
  if (argparse.HasFlag("r")) {
    package = std::move(std::unique_ptr<bridge::RnaPackage>(
        new bridge::Rnastructure("extern/rnark/data_tables/", false)));
  } else if (argparse.HasFlag("m")) {
    package = std::move(std::unique_ptr<bridge::RnaPackage>(
        new bridge::Rnark("extern/rnark/data_tables/")));
  } else {
    package = std::move(std::unique_ptr<bridge::RnaPackage>(
        new bridge::Memerna("data/")));
  }

  const auto& p = argparse.GetPositional();
  bool read_stdin = p.empty();
  std::deque<std::string> pos(p.begin(), p.end());
  if (argparse.HasFlag("e")) {
    while (1) {
      std::string seq, db;
      if (read_stdin) {
        getline(std::cin, seq);
        getline(std::cin, db);
        if (!std::cin) break;
      } else {
        if (pos.empty()) break;
        verify_expr(pos.size() >= 2, "need even number of positional args");
        seq = pos.front();
        pos.pop_front();
        db = pos.front();
        pos.pop_front();
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
        if (pos.empty()) break;
        seq = pos.front();
        pos.pop_front();
      }
      auto rna = parsing::StringToRna(seq);
      auto res = package->Fold(rna);
      printf("%d\n%s\n", res.energy, parsing::PairsToDotBracket(res.p).c_str());
    }
  }
}
