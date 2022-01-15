// Copyright 2016 Eliot Courtney.
#include <cstdio>

#include "compute/energy/load_model.h"
#include "model/context.h"

int main(int argc, char* argv[]) {
  mrna::ArgParse args(mrna::energy::ENERGY_OPTS);
  args.AddOptions(mrna::MODEL_OPTS);
  args.ParseOrExit(argc, argv);
  const auto& pos = args.GetPositional();
  verify(pos.size() == 1, "need primary sequence to fold");

  mrna::Context ctx(mrna::Primary::FromString(pos.front()),
      mrna::energy::LoadEnergyModelFromArgParse(args), ModelCfgFromArgParse(args));
  const auto res = ctx.Fold();

  printf("Energy: %d\n%s\n%s\n", res.mfe.energy, res.tb.s.ToDotBracket().c_str(),
      res.tb.ctd.ToString(res.tb.s).c_str());
}
