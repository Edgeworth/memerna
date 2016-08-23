#include <cstdio>
#include "fold/fold.h"
#include "fold/traceback.h"
#include "parsing.h"
#include "energy/energy_model.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  ArgParse argparse(energy::ENERGY_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  auto pos = argparse.GetPositional();
  verify_expr(pos.size() == 4, "need primary sequence and dp table indices");

  energy::LoadEnergyModelFromArgParse(argparse);
  fold::fold_state_t state;
  fold::Fold0(parsing::StringToRna(pos[0]), &state);

  int st = atoi(pos[1].c_str()), en = atoi(pos[2].c_str()), a = atoi(pos[3].c_str());
  fold::traceback_stack_t q;
  q.emplace(st, en, a);
  fold::TraceStructure(state.dp_table, q);

  printf("DP value at %d %d %d: %d\n  %s\n", st, en, a,
      state.dp_table[st][en][a], parsing::PairsToDotBracket(p).c_str());

  printf("Multiloop initiation params: %d %d\n", g_multiloop_hack_a, g_multiloop_hack_b);
}
