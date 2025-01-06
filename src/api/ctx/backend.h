// Copyright 2024 Eliot Courtney.
#ifndef API_CTX_BACKEND_H_
#define API_CTX_BACKEND_H_
#include <vector>

#include "api/energy/energy.h"
#include "backends/base/energy/boltz_model.h"
#include "backends/base/energy/model.h"
#include "backends/baseopt/energy/boltz_model.h"
#include "backends/baseopt/energy/model.h"
#include "backends/stack/energy/boltz_model.h"
#include "backends/stack/energy/model.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/argparse.h"

namespace mrna {

using BackendModelPtr =
    std::variant<md::base::Model::Ptr, md::base::opt::Model::Ptr, md::stack::Model::Ptr>;
using BackendBoltzModelPtr = std::variant<md::base::BoltzModel::Ptr, md::base::opt::BoltzModel::Ptr,
    md::stack::BoltzModel::Ptr>;

BackendModelPtr BackendFromArgParse(const ArgParse& args);

BackendModelPtr BackendFromBackendCfg(const BackendCfg& cfg);

constexpr BackendKind GetBackendKind(const BackendModelPtr& m) {
  auto vis = overloaded{
      [](const md::base::Model::Ptr&) { return BackendKind::BASE; },
      [](const md::base::opt::Model::Ptr&) { return BackendKind::BASEOPT; },
      [](const md::stack::Model::Ptr&) { return BackendKind::STACK; },
  };
  return std::visit(vis, m);
}

BackendModelPtr Random(BackendKind kind, uint_fast32_t seed);

// Creates the Boltzmann energy model from the given energy model.
BackendBoltzModelPtr Boltz(const BackendModelPtr& m);

// Gets energy configuration for the given energy model.
[[nodiscard]] constexpr erg::EnergyCfg BackendEnergyCfg(const BackendModelPtr& m) {
  return std::visit([&](const auto& m) { return m->cfg(); }, m);
}

void LoadPseudofreeEnergy(
    const BackendModelPtr& m, std::vector<Energy> pf_paired, std::vector<Energy> pf_unpaired);

BackendModelPtr CloneBackend(const BackendModelPtr& m);

// Returns the underlying non-Boltzmann energy model for the given Boltzmann
// energy model. Note that this may be different to the original energy model,
// i.e. m != BoltzUnderlying(Boltz(m)). For example, Boltzing an energy model turns
// off bulge loop C state calculation.
// This may create a new energy model, so it's expensive.
[[nodiscard]] inline BackendModelPtr Underlying(const BackendBoltzModelPtr& bm) {
  return std::visit([](const auto& bm) -> BackendModelPtr { return bm->m().Clone(); }, bm);
}

[[nodiscard]] inline bool CanPair(const BackendModelPtr& m, const Primary& r, int st, int en) {
  return std::visit([&](const auto& m) { return m->CanPair(r, st, en); }, m);
}

[[nodiscard]] inline erg::EnergyResult TotalEnergy(const BackendModelPtr& m, const Primary& r,
    const Secondary& s, const Ctds* given_ctd, bool build_structure = false) {
  return std::visit(
      [&](const auto& m) { return m->TotalEnergy(r, s, given_ctd, build_structure); }, m);
}

[[nodiscard]] inline erg::EnergyResult SubEnergy(const BackendModelPtr& m, const Primary& r,
    const Secondary& s, const Ctds* given_ctd, int st, int en, bool build_structure = false) {
  return std::visit(
      [&](const auto& m) { return m->SubEnergy(r, s, given_ctd, st, en, build_structure); }, m);
}

}  // namespace mrna

#endif  // API_CTX_BACKEND_H_
