// Copyright 2016 Eliot Courtney.
#include "bridge/bridge.h"

#include "bridge/memerna.h"
#include "compute/energy/load_model.h"
#include "util/macros.h"

#ifdef USE_RNASTRUCTURE
#include "bridge/rnastructure.h"
#endif  // USE_RNASTRUCTURE

namespace mrna::bridge {

std::unique_ptr<RnaPackage> RnaPackageFromArgParse(const ArgParse& args) {
  verify(args.HasFlag("r") + args.HasFlag("k") == 1, "require exactly one package flag\n%s",
      args.Usage().c_str());
  if (args.HasFlag("r")) {
#ifdef USE_RNASTRUCTURE
    return std::unique_ptr<RnaPackage>(
        new RNAstructure(args.GetOption("rnastructure-data"), false));
#else
    error("not compiled with RNAstructure");
#endif  // USE_RNASTRUCTURE

  } else {
    return std::unique_ptr<RnaPackage>(
        new Memerna(energy::LoadEnergyModelFromArgParse(args), ModelCfgFromArgParse(args)));
  }
}

}  // namespace mrna::bridge
