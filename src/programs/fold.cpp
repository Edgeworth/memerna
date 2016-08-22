#include <cstdio>
#include "parsing.h"
#include "energy/energy_model.h"
#include "fold/fold.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  ArgParse argparse(energy::ENERGY_OPTIONS);
  argparse.AddOptions(fold::FOLD_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  auto pos = argparse.GetPositional();
  verify_expr(pos.size() == 1, "need primary sequence to fold");

  energy::LoadEnergyModelFromArgParse(argparse);
  auto fold_fn = fold::FoldFunctionFromArgParse(argparse);
  auto frna = fold_fn(parsing::StringToRna(pos.front()), nullptr);

  printf("Energy: %d\n%s\n", frna.energy, parsing::PairsToDotBracket(frna.p).c_str());
}
