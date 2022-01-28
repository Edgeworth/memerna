// Copyright 2016 E.
#include "bridge/bridge.h"

#include "bridge/memerna.h"
#include "compute/energy/energy.h"
#include "compute/energy/model.h"
#include "context/config.h"
#include "util/error.h"

#ifdef USE_RNASTRUCTURE
#include "bridge/rnastructure.h"
#endif  // USE_RNASTRUCTURE

namespace mrna::bridge {

void RegisterOpts(ArgParse* args) {
  ::mrna::RegisterOpts(args);
  args->RegisterOpt(OPT_USE_RNASTRUCTURE);
  args->RegisterOpt(OPT_RNASTRUCTURE_DATA);
  args->RegisterOpt(OPT_USE_MEMERNA);
}

std::unique_ptr<RnaPackage> RnaPackage::FromArgParse(const ArgParse& args) {
  verify(args.Has(OPT_USE_RNASTRUCTURE) + args.Has(OPT_USE_MEMERNA) == 1,
      "require exactly one package flag\n%s", args.Usage().c_str());
  if (args.Has(OPT_USE_RNASTRUCTURE)) {
#ifdef USE_RNASTRUCTURE
    return std::unique_ptr<RnaPackage>(new RNAstructure(args.Get(OPT_RNASTRUCTURE_DATA), false));
#else
    error("not compiled with RNAstructure");
#endif  // USE_RNASTRUCTURE

  } else {
    return std::unique_ptr<RnaPackage>(new Memerna(Memerna::FromArgParse(args)));
  }
}

}  // namespace mrna::bridge
