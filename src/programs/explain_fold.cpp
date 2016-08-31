#include <cstdio>
#include "fold/fold.h"
#include "fold/traceback.h"
#include "parsing.h"
#include "energy/energy_model.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  ArgParse argparse(energy::ENERGY_OPTIONS);
  argparse.AddOptions(fold::FOLD_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  auto pos = argparse.GetPositional();
  verify_expr(pos.size() == 4, "need primary sequence and dp table indices");

  energy::LoadEnergyModelFromArgParse(argparse);
  auto r = parsing::StringToPrimary(pos[0]);
  fold::fold_state_t state;
  auto fold_fn = fold::FoldFunctionFromArgParse(argparse);
  fold_fn(r, &state);

  int st = atoi(pos[1].c_str()), en = atoi(pos[2].c_str()), a = atoi(pos[3].c_str());
  fold::traceback_stack_t q;
  q.emplace(st, en, a);
  auto computed = fold::TraceStructure(r, state.dp_table, state.ext_table, q);

  printf("DP value at %d %d %d: %d\n  %s\n", st, en, a,
      state.dp_table[st][en][a], parsing::PairsToDotBracket(computed.p).c_str());
}
