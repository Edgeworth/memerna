// Copyright 2026 Eliot Courtney.
#include "backends/canon/energy/model.h"

#include <string>

namespace mrna::md::canon {

bool Model::IsSupported(const BackendCfg& backend_cfg, const erg::EnergyCfg& /*cfg*/,
    const erg::PseudofreeCfg& /*pf*/, std::string* reason) {
  if (backend_cfg.energy_model == erg::EnergyModelKind::T22) {
    if (reason) *reason = "T22 energy model not supported by canon backend";
    return false;
  }
  return true;
}

}  // namespace mrna::md::canon
