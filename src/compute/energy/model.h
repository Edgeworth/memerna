// Copyright 2022 Eliot Courtney.
#ifndef COMPUTE_ENERGY_MODEL_H_
#define COMPUTE_ENERGY_MODEL_H_
#include <cstdint>
#include <iosfwd>
#include <variant>

#include "compute/energy/energy.h"
#include "compute/energy/energy_cfg.h"
#include "compute/energy/t04/boltz_model.h"
#include "compute/energy/t04/model.h"
#include "compute/energy/t22/boltz_model.h"
#include "compute/energy/t22/model.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/argparse.h"
#include "util/error.h"
#include "util/util.h"

namespace mrna::erg {

using EnergyModelPtr = std::variant<t04::ModelPtr, t22::ModelPtr>;
using BoltzEnergyModelPtr = std::variant<t04::BoltzModelPtr, t22::BoltzModelPtr>;

enum class ModelKind {
  T04_LIKE,
  T22_LIKE,
};

std::istream& operator>>(std::istream& is, ModelKind& kind);
std::ostream& operator<<(std::ostream& os, ModelKind kind);

inline EnergyModelPtr FromArgParse(const ArgParse& args) {
  auto kind = args.Get<ModelKind>(OPT_ENERGY_MODEL);
  switch (kind) {
  case ModelKind::T04_LIKE: return t04::Model::FromArgParse(args);
  case ModelKind::T22_LIKE: return t22::Model::FromArgParse(args);
  default: bug();
  }
}

inline EnergyModelPtr Random(ModelKind kind, uint_fast32_t seed) {
  switch (kind) {
  case ModelKind::T04_LIKE: return t04::Model::Random(seed);
  case ModelKind::T22_LIKE: return t22::Model::Random(seed);
  default: bug();
  }
}

// TODO(0): can simplify all these visits via ModelMixin.

// Creates the Boltzmann energy model from the given energy model.
inline BoltzEnergyModelPtr Boltz(const EnergyModelPtr& em) {
  auto vis = overloaded{
      [](const t04::ModelPtr& em) -> BoltzEnergyModelPtr { return t04::BoltzModel::Create(em); },
      [](const t22::ModelPtr& em) -> BoltzEnergyModelPtr { return t22::BoltzModel::Create(em); },
  };
  return std::visit(vis, em);
}

// Returns the underlying non-Boltzmann energy model for the given Boltzmann
// energy model. Note that this may be different to the original energy model,
// i.e. em != BoltzUnderlying(Boltz(em)). For example, Boltzing an energy model turns
// off bulge loop C state calculation.
// This may create a new energy model, so it's expensive.
inline EnergyModelPtr Underlying(const BoltzEnergyModelPtr& bem) {
  return std::visit([](const auto& bem) -> EnergyModelPtr { return bem->em().Clone(); }, bem);
}

inline ModelKind Kind(const EnergyModelPtr& em) {
  auto vis = overloaded{
      [](const t04::ModelPtr&) { return ModelKind::T04_LIKE; },
      [](const t22::ModelPtr&) { return ModelKind::T22_LIKE; },
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

#endif  // COMPUTE_ENERGY_MODEL_H_
