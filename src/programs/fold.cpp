#include <cstdio>
#include <climits>
#include "fold/suboptimal0.h"
#include "parsing.h"
#include "energy/energy_model.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  ArgParse argparse(energy::ENERGY_OPTIONS);
  argparse.AddOptions(fold::FOLD_OPTIONS);
  argparse.AddOptions({
      {"subopt-delta", ArgParse::option_t("maximum energy delta").Arg("-1")},
      {"subopt-max", ArgParse::option_t("maximum number of reported structures").Arg("-1")}
  });
  argparse.ParseOrExit(argc, argv);
  auto pos = argparse.GetPositional();
  verify_expr(pos.size() == 1, "need primary sequence to fold");

  energy::LoadEnergyModelFromArgParse(argparse);
  auto fold_fn = fold::FoldFunctionFromArgParse(argparse);
  fold::fold_state_t state;
  auto r = parsing::StringToPrimary(pos.front());
  auto computed = fold_fn(r, &state);

  printf("Energy: %d\n%s\n", computed.energy, parsing::PairsToDotBracket(computed.p).c_str());

  energy_t energy_delta = atoi(argparse.GetOption("subopt-delta").c_str());
  int max_structures = atoi(argparse.GetOption("subopt-max").c_str());
  if (energy_delta >= 0 || max_structures > 0) {
    if (energy_delta < 0) energy_delta = constants::CAP_E;
    if (max_structures < 0) max_structures = INT_MAX;
    auto structures = fold::SuboptimalTraceback0(
        computed.energy + energy_delta, max_structures, r, state.dp_table, state.ext_table);
    printf("%zu suboptimal structures:\n", structures.size());
    for (const auto& subopt : structures) {
      printf("%d %s\n", subopt.energy, parsing::PairsToDotBracket(subopt.p).c_str());
    }
  }
}
