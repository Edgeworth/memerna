// Copyright 2016 E.
#include "bridge/bridge.h"

#include "bridge/memerna.h"
#include "bridge/rnastructure.h"
#include "energy/load_model.h"

namespace mrna {
namespace bridge {

std::unique_ptr<RnaPackage> RnaPackageFromArgParse(const ArgParse& argparse) {
  verify(argparse.HasFlag("r") + argparse.HasFlag("k") == 1, "require exactly one package flag\n%s",
      argparse.Usage().c_str());
  if (argparse.HasFlag("r")) {
    return std::unique_ptr<RnaPackage>(new RNAstructure(argparse.GetOption("data-path"), false));
  } else {
    return std::unique_ptr<RnaPackage>(new Memerna(
        energy::LoadEnergyModelFromArgParse(argparse), ContextOptionsFromArgParse(argparse)));
  }
}
}  // namespace bridge
}  // namespace mrna
