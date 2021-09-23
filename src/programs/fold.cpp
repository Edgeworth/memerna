// Copyright 2016 Eliot Courtney.
#include <cstdio>

#include "context.h"
#include "energy/load_model.h"
#include "parsing.h"

int main(int argc, char* argv[]) {
  mrna::ArgParse args(mrna::energy::ENERGY_OPTIONS);
  args.AddOptions(mrna::CONTEXT_OPTIONS);
  args.ParseOrExit(argc, argv);
  const auto& pos = args.GetPositional();
  verify(pos.size() == 1, "need primary sequence to fold");

  mrna::Context ctx(mrna::parsing::StringToPrimary(pos.front()),
      mrna::energy::LoadEnergyModelFromArgParse(args), ContextOptionsFromArgParse(args));
  const auto computed = ctx.Fold();

  printf("Energy: %d\n%s\n%s\n", computed.energy,
      mrna::parsing::PairsToDotBracket(computed.s.p).c_str(),
      mrna::parsing::ComputedToCtdString(computed).c_str());
}
