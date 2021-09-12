// Copyright 2016 Eliot Courtney.
#include <cstdio>

#include "context.h"
#include "energy/load_model.h"
#include "parsing.h"

using namespace mrna;

int main(int argc, char* argv[]) {
  ArgParse argparse(energy::ENERGY_OPTIONS);
  argparse.AddOptions(CONTEXT_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  const auto& pos = argparse.GetPositional();
  verify(pos.size() == 1, "need primary sequence to fold");

  Context ctx(parsing::StringToPrimary(pos.front()), energy::LoadEnergyModelFromArgParse(argparse),
      ContextOptionsFromArgParse(argparse));
  const auto computed = ctx.Fold();

  printf("Energy: %d\n%s\n%s\n", computed.energy, parsing::PairsToDotBracket(computed.s.p).c_str(),
      parsing::ComputedToCtdString(computed).c_str());
}
