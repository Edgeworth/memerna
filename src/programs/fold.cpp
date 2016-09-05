#include <cstdio>
#include "fold/suboptimal0.h"
#include "parsing.h"
#include "energy/load_model.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  ArgParse argparse(energy::ENERGY_OPTIONS);
  argparse.AddOptions(fold::FOLD_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  const auto pos = argparse.GetPositional();
  verify_expr(pos.size() == 1, "need primary sequence to fold");

  fold::Context ctx(
      parsing::StringToPrimary(pos.front()),
      energy::LoadEnergyModelFromArgParse(argparse),
      fold::ContextOptionsFromArgParse(argparse));
  const auto computed = ctx.Fold();

  printf("Energy: %d\n%s\n%s\n",
      computed.energy, parsing::PairsToDotBracket(computed.s.p).c_str(),
      parsing::ComputedToCtdString(computed).c_str());

  if (ctx.GetOptions().subopt_energy >= 0 || ctx.GetOptions().subopt_num > 0) {
    const auto structures = ctx.Suboptimal();
    printf("%zu suboptimal structures:\n", structures.size());
    for (const auto& subopt : structures) {
      printf("%d %s\n", subopt.energy, parsing::ComputedToCtdString(subopt).c_str());
    }
  }
}
