// Copyright 2016 E.
#include "bridge/bridge.h"

#include "bridge/memerna.h"
#include "bridge/rnastructure.h"
#include "energy/load_model.h"

namespace mrna {
namespace bridge {

std::unique_ptr<RnaPackage> RnaPackageFromArgParse(const ArgParse& args) {
  verify(args.HasFlag("r") + args.HasFlag("k") == 1, "require exactly one package flag\n%s",
      args.Usage().c_str());
  if (args.HasFlag("r")) {
    return std::unique_ptr<RnaPackage>(
        new RNAstructure(args.GetOption("rnastructure-data"), false));
  } else {
    return std::unique_ptr<RnaPackage>(
        new Memerna(energy::LoadEnergyModelFromArgParse(args), ContextOptionsFromArgParse(args)));
  }
}

}  // namespace bridge
}  // namespace mrna
