// Copyright 2016 E.
#include "bridge/bridge.h"

#include "bridge/memerna.h"
#include "compute/energy/model.h"
#include "model/config.h"
#include "util/error.h"

#ifdef USE_RNASTRUCTURE
#include "bridge/rnastructure.h"
#endif  // USE_RNASTRUCTURE

namespace mrna::bridge {

std::unique_ptr<RnaPackage> RnaPackage::FromArgParse(const ArgParse& args) {
  verify(args.HasFlag("r") + args.HasFlag("m") == 1, "require exactly one package flag\n%s",
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
        new Memerna(energy::EnergyModel::FromArgParse(args), ModelCfg::FromArgParse(args)));
  }
}

}  // namespace mrna::bridge
