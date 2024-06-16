// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_STACK_ENERGY_MODEL_H_
#define BACKENDS_STACK_ENERGY_MODEL_H_

#include <fmt/core.h>

#include <cassert>
#include <cmath>
#include <deque>
#include <memory>
#include <random>
#include <string>

#include "api/energy/energy.h"
#include "backends/common/base/model_base.h"
#include "backends/common/base/parse.h"
#include "backends/common/model_mixin.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "model/structure.h"

namespace mrna::md::stack {

using mrna::erg::EnergyResult;

class Model : public base::ModelBase, public ModelMixin<Model> {
 public:
  static constexpr auto KIND = BackendKind::STACK;

  Energy penultimate_stack[4][4][4][4] = {};

  // Pseudofree energies. Ignored if empty.
  std::vector<Energy> pf_paired;
  std::vector<Energy> pf_unpaired;
  // Cumulative sum of size N+1 (first element is nothing).
  std::vector<Energy> pf_unpaired_cum;

  Energy Hairpin(const Primary& r, int st, int en, std::unique_ptr<Structure>* s = nullptr) const;
  Energy Bulge(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const;
  Energy InternalLoop(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const;
  Energy TwoLoop(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const;
  Energy MultiloopEnergy(const Primary& r, const Secondary& s, int st, int en,
      std::deque<int>* branches, bool use_given_ctds, Ctds* ctd,
      std::unique_ptr<Structure>* sstruc = nullptr) const;

  [[nodiscard]] constexpr Energy PfUnpaired(int n) const {
    if (pf_unpaired.empty()) return ZERO_E;
    return pf_unpaired[n];
  }

  // Inclusive range, unlike pf_unpaired_cum directly.
  [[nodiscard]] constexpr Energy PfUnpairedCum(int st, int en) const {
    if (pf_unpaired.empty()) return ZERO_E;
    return pf_unpaired_cum[en + 1] - pf_unpaired_cum[st];
  }

  [[nodiscard]] constexpr Energy PfPaired(int st, int en) const {
    if (pf_paired.empty()) return ZERO_E;
    return pf_paired[st] + pf_paired[en];
  }

  // Computes the penalty for a stack of the given length, ending at (ist, ien).
  // Handles bulge loops.
  [[nodiscard]] constexpr Energy StackPenalty(const Primary& r, const Secondary& s, int ost,
      int oen, int ist, int ien, std::unique_ptr<Structure>* struc = nullptr) const;

  // ModelMixin:
  EnergyResult SubEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd, int st,
      int en, bool build_structure = false) const;
  EnergyResult TotalEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd,
      bool build_structure = false) const;

  bool IsValid(std::string* reason = nullptr) const { return base::ModelIsValid(*this, reason); }

  void VerifyValidFor(const Primary& r) const {
    if (!pf_paired.empty())
      verify(pf_paired.size() == r.size(), "pseudofree paired must be same length as seq");
    if (!pf_unpaired.empty())
      verify(pf_unpaired.size() == r.size(), "pseudofree unpaired must be same length as seq");
  }

  void LoadPseudofreeEnergy(std::vector<Energy> paired, std::vector<Energy> unpaired);

  void LoadFromModelPath(const std::string& path);
  void LoadRandom(std::mt19937& eng);

 private:
  friend class ModelMixin<Model>;

  void VerifyLengths(const Primary& r, const Secondary& s, const Ctds* given_ctd) const;

  // This is private to prevent construction on the stack, since this structure is large.
  Model() = default;

  Energy SubEnergyInternal(const Primary& r, const Secondary& s, int st, int en, int stack_st,
      int stack_en, bool use_given_ctds, Ctds* ctd, std::unique_ptr<Structure>* struc) const;
};

}  // namespace mrna::md::stack

#endif  // BACKENDS_STACK_ENERGY_MODEL_H_
