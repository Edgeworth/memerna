#include "energy/load_model.h"
#include "energy/energy.h"
#include "fold/context.h"
#include "bridge/bridge.h"
#include "bridge/rnastructure.h"
#include "bridge/memerna.h"

namespace memerna {
namespace bridge {

std::unique_ptr<RnaPackage> RnaPackageFromArgParse(const ArgParse& argparse) {
  verify_expr(
      argparse.HasFlag("r") + argparse.HasFlag("k") == 1,
      "require exactly one package flag\n%s", argparse.Usage().c_str());
  if (argparse.HasFlag("r")) {
    return std::unique_ptr<RnaPackage>(new Rnastructure("extern/miles_rnastructure/data_tables/", false));
  } else {
    return std::unique_ptr<RnaPackage>(new Memerna(
        energy::LoadEnergyModelFromArgParse(argparse), fold::ContextOptionsFromArgParse(argparse)));
  }
}

}
}
