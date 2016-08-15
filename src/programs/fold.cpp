#include <cstdio>
#include "argparse.h"
#include "base.h"
#include "parsing.h"
#include "fold/fold.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  ArgParse argparse({
      {"print-interval", ArgParse::option_t("status update every n seconds").Arg("60")}
  });
  auto ret = argparse.Parse(argc, argv);
  verify_expr(
      ret.size() == 0,
      "%s\n%s\n", ret.c_str(), argparse.Usage().c_str());
  LoadEnergyModelFromDataDir("data");
  SetRna(parsing::StringToRna(argv[1]));
  energy_t e = fold::Fold1();
  printf("Energy: %d\n%s\n", e, parsing::PairsToDotBracket(p).c_str());
}
