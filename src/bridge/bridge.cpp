// Copyright 2016 Eliot Courtney.
#include "bridge/bridge.h"

#include "bridge/memerna.h"
#include "bridge/rnastructure.h"
#include "energy/load_model.h"

namespace memerna {
namespace bridge {

std::unique_ptr<RnaPackage> RnaPackageFromArgParse(const ArgParse& argparse) {
  verify_expr(argparse.HasFlag("r") + argparse.HasFlag("k") == 1,
      "require exactly one package flag\n%s", argparse.Usage().c_str());
  if (argparse.HasFlag("r")) {
    return std::unique_ptr<RnaPackage>(
        new Rnastructure("extern/miles_rnastructure/data_tables/", false));
  } else {
    return std::unique_ptr<RnaPackage>(new Memerna(
        energy::LoadEnergyModelFromArgParse(argparse), ContextOptionsFromArgParse(argparse)));
  }
}
}  // namespace bridge
}  // namespace memerna
