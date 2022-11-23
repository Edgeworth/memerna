#ifndef COMPUTE_ENERGY_MODEL_H_
#define COMPUTE_ENERGY_MODEL_H_
#include <memory>
#include <variant>

#include "compute/energy/t04/boltz_model.h"
#include "compute/energy/t04/model.h"
#include "model/constants.h"
#include "model/ctd.h"

namespace mrna::energy {

using EnergyModelPtr = std::variant<t04::ModelPtr>;
using BoltzEnergyModelPtr = std::variant<t04::BoltzModelPtr>;

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
  case ModelKind::T22_LIKE:
    // TODO(0): implement.
    bug();
  default: bug();
  }
}

inline EnergyModelPtr Random(ModelKind kind, uint_fast32_t seed) {
  switch (kind) {
  case ModelKind::T04_LIKE: return t04::Model::Random(seed);
  case ModelKind::T22_LIKE:
    // TODO(0): implement.
    bug();
  default: bug();
  }
}

inline BoltzEnergyModelPtr Boltz(const EnergyModelPtr& em) {
  return std::visit([](const t04::ModelPtr& em) { return t04::BoltzModel::Create(em); }, em);
}

inline EnergyModelPtr Unboltz(const BoltzEnergyModelPtr& bem) {
  return std::visit([](const t04::BoltzModelPtr& bem) { return bem->em(); }, em);
}

inline ModelKind Kind(const EnergyModelPtr& model) {
  return std::visit([](const t04::ModelPtr&) { return ModelKind::T04_LIKE; }, model);
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

}  // namespace mrna::energy

#endif  // COMPUTE_ENERGY_MODEL_H_
