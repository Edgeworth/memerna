#include <cstdio>
#include <bridge/bridge.h>
#include "parsing.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  ArgParse argparse;
  argparse.AddOptions(MEMERNA_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  auto pos = argparse.GetPositional();
  verify_expr(pos.size() == 1, "need primary sequence to fold");

  LoadEnergyModelFromDataDir("data/");
  auto fold_fn = FoldFunctionFromArgParse(argparse);
  auto frna = fold_fn(parsing::StringToRna(pos.front()));

  printf("Energy: %d\n%s\n", frna.energy, parsing::PairsToDotBracket(frna.p).c_str());
}
