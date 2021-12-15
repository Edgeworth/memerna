// Copyright 2016 Eliot Courtney.
#include <cstdio>

#include "energy/load_model.h"
#include "model/context.h"
#include "model/parsing.h"

int main(int argc, char* argv[]) {
  mrna::ArgParse args(mrna::energy::ENERGY_OPTIONS);
  args.AddOptions(mrna::CONTEXT_OPTIONS);
  args.ParseOrExit(argc, argv);
  const auto& pos = args.GetPositional();
  verify(pos.size() == 1, "need primary sequence to fold");

  mrna::Context ctx(mrna::StringToPrimary(pos.front()),
      mrna::energy::LoadEnergyModelFromArgParse(args), ContextOptionsFromArgParse(args));
  const auto computed = ctx.Fold();

  printf("Energy: %d\n%s\n%s\n", computed.energy, mrna::PairsToDotBracket(computed.s.p).c_str(),
      mrna::ComputedToCtdString(computed).c_str());
}
