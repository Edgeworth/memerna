// Copyright 2016 E.
#include <cstdio>

#include "context.h"
#include "energy/load_model.h"
#include "parsing.h"

int main(int argc, char* argv[]) {
  mrna::ArgParse argparse(mrna::energy::ENERGY_OPTIONS);
  argparse.AddOptions(mrna::CONTEXT_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  const auto& pos = argparse.GetPositional();
  verify(pos.size() == 1, "need primary sequence to fold");

  mrna::Context ctx(mrna::parsing::StringToPrimary(pos.front()),
      mrna::energy::LoadEnergyModelFromArgParse(argparse), ContextOptionsFromArgParse(argparse));
  const auto computed = ctx.Fold();

  printf("Energy: %d\n%s\n%s\n", computed.energy,
      mrna::parsing::PairsToDotBracket(computed.s.p).c_str(),
      mrna::parsing::ComputedToCtdString(computed).c_str());
}
