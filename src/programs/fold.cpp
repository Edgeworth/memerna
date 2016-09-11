#include <cstdio>
#include "fold/context.h"
#include "energy/load_model.h"
#include "parsing.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  ArgParse argparse(energy::ENERGY_OPTIONS);
  argparse.AddOptions(fold::FOLD_OPTIONS);
  argparse.AddOptions({
      {"subopt-delta", ArgParse::option_t("maximum energy delta from minimum").Arg("-1")},
      {"subopt-num", ArgParse::option_t("maximum number of reported structures").Arg("-1")}
  });
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

  energy_t subopt_delta = atoi(argparse.GetOption("subopt-delta").c_str());
  int subopt_num = atoi(argparse.GetOption("subopt-num").c_str());
  if (subopt_delta >= 0 || subopt_num > 0) {
    const auto structures = ctx.Suboptimal(subopt_delta, subopt_num);
    printf("%zu suboptimal structures:\n", structures.size());
    for (const auto& subopt : structures) {
      printf("%d %s\n", subopt.energy, parsing::ComputedToCtdString(subopt).c_str());
    }
  }
}
