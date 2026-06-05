// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_CANON_MFE_MFE_OPT_H_
#define BACKENDS_CANON_MFE_MFE_OPT_H_

#include <string>

#include "backends/canon/dp.h"
#include "backends/canon/energy/model.h"
#include "model/primary.h"

namespace mrna::md::canon {

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
class MfeOpt {
 public:
  static bool IsSupported(const PseudofreeCfgT& pf, std::string* reason = nullptr);

  static void Run(const Primary& r, const Model::Ptr& m, DpState& state, const PseudofreeCfgT& pf);
};

}  // namespace mrna::md::canon

#endif  // BACKENDS_CANON_MFE_MFE_OPT_H_
