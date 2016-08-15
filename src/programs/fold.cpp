#include <cstdio>
#include "base.h"
#include "parsing.h"
#include "energy/structure.h"
#include "fold/fold.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  verify_expr(argc == 2, "requires one argument");
  LoadEnergyModelFromDataDir("data");
  SetRna(parsing::StringToRna(argv[1]));
  energy_t e = fold::Fold1();
  printf("Energy: %d\n%s\n", e, parsing::PairsToDotBracket(p).c_str());
}
