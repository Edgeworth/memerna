#include <iostream>
#include <memory>
#include "parsing.h"
#include "argparse.h"
#include "bridge/bridge.h"

// efn, fold for the three progs + brute fold?

int main(int argc, char* argv[]) {
  memerna::ArgParse argparse({
      {"v", "be verbose (if possible)"},
      {"r", "rnastructure"},
      {"m", "rnark"},
      {"k", "memerna"},
      {"e", "run efn"},
      {"f", "run fold"}
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

  std::unique_ptr<memerna::bridge::RnaPackage> package;
  if (argparse.HasFlag("r")) {
    package = std::move(std::unique_ptr<memerna::bridge::RnaPackage>(
        new memerna::bridge::Rnastructure("extern/rnark/data_tables")));
  } else if (argparse.HasFlag("m")) {
    package = std::move(std::unique_ptr<memerna::bridge::RnaPackage>(
        new memerna::bridge::Rnark("extern/rnark/data_tables")));
  } else {
    package = std::move(std::unique_ptr<memerna::bridge::RnaPackage>(
        new memerna::bridge::Memerna("data")));
  }

  if (argparse.HasFlag("e")) {
    while (1) {
      std::string seq, db;
      getline(std::cin, seq);
      getline(std::cin, db);
      if (!std::cin) break;
      auto frna = memerna::parsing::ParseDotBracketRna(seq, db);
      auto res = package->Efn(frna, false);
      printf("%d\n%s", res.energy, res.desc.c_str());
    }
  } else {
    while (1) {
      std::string seq;
      getline(std::cin, seq);
      if (!std::cin) break;
      auto rna = memerna::parsing::StringToRna(seq);
      auto res = package->Fold(rna, false);
      printf("%d\n%s\n%s", res.energy,
          memerna::parsing::PairsToDotBracket(res.frna.p).c_str(), res.desc.c_str());
    }
  }
}
