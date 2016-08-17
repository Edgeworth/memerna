#include <cstdio>
#include <bridge/bridge.h>
#include "argparse.h"
#include "base.h"
#include "parsing.h"
#include "fold/fold.h"

using namespace memerna;

int main(int argc, char* argv[]) {
  ArgParse argparse({
      {"print-interval", ArgParse::option_t("status update every n seconds").Arg("60")},
  });
  argparse.AddOptions(bridge::MEMERNA_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  auto pos = argparse.GetPositional();
  verify_expr(pos.size() == 1, "need primary sequence to fold");

  auto memerna = bridge::MemernaFromArgParse(argparse);
  auto frna = memerna->Fold(parsing::StringToRna(pos.front()));

  printf("Energy: %d\n%s\n", frna.energy, parsing::PairsToDotBracket(frna.p).c_str());
}
