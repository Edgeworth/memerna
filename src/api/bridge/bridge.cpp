// Copyright 2016 Eliot Courtney.
#include "api/bridge/bridge.h"

#include "api/bridge/memerna.h"
#include "api/ctx/ctx_cfg.h"
#include "util/error.h"

#ifdef USE_RNASTRUCTURE
#include "api/bridge/rnastructure.h"
#endif  // USE_RNASTRUCTURE

namespace mrna::bridge {

void RegisterOpts(ArgParse* args) {
  ::mrna::RegisterOpts(args);
  args->RegisterOpt(OPT_USE_RNASTRUCTURE);
  args->RegisterOpt(OPT_RNASTRUCTURE_DATA);
  args->RegisterOpt(OPT_USE_MEMERNA);
}

std::unique_ptr<RnaPackage> RnaPackage::FromArgParse(const ArgParse& args) {
  verify(args.Get<bool>(OPT_USE_RNASTRUCTURE) + args.Get<bool>(OPT_USE_MEMERNA) == 1,
      "require exactly one package flag\n{}", args.Usage());
  if (args.Get<bool>(OPT_USE_RNASTRUCTURE)) {
#ifdef USE_RNASTRUCTURE
    return std::unique_ptr<RnaPackage>(new RNAstructure(RNAstructure::FromArgParse(args)));
#else
    error("not compiled with RNAstructure");
#endif  // USE_RNASTRUCTURE

  } else {
    return std::unique_ptr<RnaPackage>(new Memerna(Memerna::FromArgParse(args)));
  }
}

}  // namespace mrna::bridge
