// Copyright 2022 Eliot Courtney.
#ifndef API_ENERGY_MODEL_H_
#define API_ENERGY_MODEL_H_

#include <cstdint>
#include <iosfwd>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "api/energy/energy.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "models/t04/energy/boltz_model.h"
#include "models/t04/energy/model.h"
#include "models/t22/energy/boltz_model.h"
#include "models/t22/energy/model.h"
#include "util/argparse.h"
#include "util/error.h"
#include "util/string.h"
#include "util/util.h"

namespace mrna::erg {

using EnergyModelPtr = std::variant<md::t04::Model::Ptr, md::t22::Model::Ptr>;
using BoltzEnergyModelPtr = std::variant<md::t04::BoltzModel::Ptr, md::t22::BoltzModel::Ptr>;

enum class ModelKind {
  T04_LIKE,
  T22_LIKE,
};

std::istream& operator>>(std::istream& is, ModelKind& kind);
std::ostream& operator<<(std::ostream& os, ModelKind kind);

inline EnergyModelPtr FromArgParse(const ArgParse& args) {
  auto kind = args.Get<ModelKind>(OPT_ENERGY_MODEL);
  switch (kind) {
  case ModelKind::T04_LIKE: return md::t04::Model::FromArgParse(args);
  case ModelKind::T22_LIKE: return md::t22::Model::FromArgParse(args);
  }
  unreachable();
}

inline EnergyModelPtr FromDir(const std::string& data_dir, const std::string& model) {
  auto kind = Conv<ModelKind>(model);
  auto path = ModelPath(data_dir, model);
  switch (kind) {
  case ModelKind::T04_LIKE: return md::t04::Model::FromModelPath(path);
  case ModelKind::T22_LIKE: return md::t22::Model::FromModelPath(path);
  }
  unreachable();
}

inline EnergyModelPtr Random(ModelKind kind, uint_fast32_t seed) {
  switch (kind) {
  case ModelKind::T04_LIKE: return md::t04::Model::Random(seed);
  case ModelKind::T22_LIKE: return md::t22::Model::Random(seed);
  }
  unreachable();
}

// Creates the Boltzmann energy model from the given energy model.
inline BoltzEnergyModelPtr Boltz(const EnergyModelPtr& em) {
  auto vis = overloaded{
      [](const md::t04::Model::Ptr& em) -> BoltzEnergyModelPtr {
        return md::t04::BoltzModel::Create(em);
      },
      [](const md::t22::Model::Ptr& em) -> BoltzEnergyModelPtr {
        return md::t22::BoltzModel::Create(em);
      },
  };
  return std::visit(vis, em);
}

// Gets energy configuration for the given energy model.
inline constexpr EnergyCfg ModelEnergyCfg(const EnergyModelPtr& em) {
  return std::visit([&](const auto& em) { return em->cfg; }, em);
}

// Creates the Boltzmann energy model from the given energy model.
inline void LoadPseudofreeEnergy(
    const EnergyModelPtr& em, std::vector<Energy> pf_paired, std::vector<Energy> pf_unpaired) {
  auto vis = overloaded{
      [&](const auto& em) mutable {
        return em->LoadPseudofreeEnergy(std::move(pf_paired), std::move(pf_unpaired));
      }
  };
  std::visit(vis, em);
}

// Returns the underlying non-Boltzmann energy model for the given Boltzmann
// energy model. Note that this may be different to the original energy model,
// i.e. em != BoltzUnderlying(Boltz(em)). For example, Boltzing an energy model turns
// off bulge loop C state calculation.
// This may create a new energy model, so it's expensive.
inline EnergyModelPtr Underlying(const BoltzEnergyModelPtr& bem) {
  return std::visit([](const auto& bem) -> EnergyModelPtr { return bem->em().Clone(); }, bem);
}

inline constexpr ModelKind Kind(const EnergyModelPtr& em) {
  auto vis = overloaded{
      [](const md::t04::Model::Ptr&) { return ModelKind::T04_LIKE; },
      [](const md::t22::Model::Ptr&) { return ModelKind::T22_LIKE; },
  };
  return std::visit(vis, em);
}

inline bool CanPair(const EnergyModelPtr& em, const Primary& r, int st, int en) {
  return std::visit([&](const auto& em) { return em->CanPair(r, st, en); }, em);
}

inline EnergyResult TotalEnergy(const EnergyModelPtr& em, const Primary& r, const Secondary& s,
    const Ctds* given_ctd, bool build_structure = false) {
  return std::visit(
      [&](const auto& em) { return em->TotalEnergy(r, s, given_ctd, build_structure); }, em);
}

inline EnergyResult SubEnergy(const EnergyModelPtr& em, const Primary& r, const Secondary& s,
    const Ctds* given_ctd, int st, int en, bool build_structure = false) {
  return std::visit(
      [&](const auto& em) { return em->SubEnergy(r, s, given_ctd, st, en, build_structure); }, em);
}

}  // namespace mrna::erg

#endif  // API_ENERGY_MODEL_H_
