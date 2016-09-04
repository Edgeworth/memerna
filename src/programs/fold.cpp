#include <cstdio>
#include "fold/suboptimal0.h"
#include "parsing.h"
#include "energy/load_model.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  ArgParse argparse(energy::ENERGY_OPTIONS);
  argparse.AddOptions(fold::FOLD_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  auto pos = argparse.GetPositional();
  verify_expr(pos.size() == 1, "need primary sequence to fold");

  fold::Context ctx(
      parsing::StringToPrimary(pos.front()),
      energy::LoadEnergyModelFromArgParse(argparse),
      fold::ContextOptionsFromArgParse(argparse));
  auto computed = ctx.Fold();

  printf("Energy: %d\n%s\n%s\n",
      computed.energy, parsing::PairsToDotBracket(computed.s.p).c_str(),
      parsing::ComputedToCtdString(computed).c_str());

  energy_t energy_delta = atoi(argparse.GetOption("subopt-delta").c_str());
  int max_structures = atoi(argparse.GetOption("subopt-max").c_str());
  if (energy_delta >= 0 || max_structures > 0) {
    auto structures = ctx.Suboptimal();
    printf("%zu suboptimal structures:\n", structures.size());
    for (const auto& subopt : structures) {
      printf("%d %s\n", subopt.energy, parsing::ComputedToCtdString(subopt).c_str());
    }
  }
}
