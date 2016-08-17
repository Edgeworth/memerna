#include <cstdio>
#include "argparse.h"
#include "base.h"
#include "parsing.h"
#include "fold/fold.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  ArgParse argparse({
                        {"print-interval", ArgParse::option_t("status update every n seconds").Arg("60")},
                        {"alg",            ArgParse::option_t("which algorithm (slow, 1, 2)").Arg("slow")}
                    });
  auto ret = argparse.Parse(argc, argv);
  verify_expr(
      ret.size() == 0,
      "%s\n%s\n", ret.c_str(), argparse.Usage().c_str());
  auto pos = argparse.GetPositional();
  verify_expr(pos.size() == 1, "need primary sequence to fold");
  LoadEnergyModelFromDataDir("data");
  SetRna(parsing::StringToRna(pos.front()));
  energy_t e;
  auto opt = argparse.GetOption("alg");
  if (opt == "slow")
    e = fold::FoldSlow();
  else if (opt == "1")
    e = fold::Fold1();
  else if (opt == "2")
    e = fold::Fold2();
  printf("Energy: %d\n%s\n", e, parsing::PairsToDotBracket(p).c_str());
}
