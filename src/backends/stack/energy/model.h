// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_STACK_ENERGY_MODEL_H_
#define BACKENDS_STACK_ENERGY_MODEL_H_

#include <fmt/core.h>

#include <cassert>
#include <deque>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "api/energy/energy.h"
#include "backends/common/base/model_base.h"
#include "backends/common/base/parse.h"
#include "backends/common/base/pseudofree_model.h"
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

  base::PseudofreeModel pf;
  Energy penultimate_stack[4][4][4][4] = {};

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

  // Computes the penalty for a stack of the given length, ending at (ist, ien).
  // Handles bulge loops.
  [[nodiscard]] constexpr Energy StackPenalty(const Primary& r, const Secondary& s, int ost,
      int oen, int ist, int ien, std::unique_ptr<Structure>* struc = nullptr) const;

  // ModelMixin:
  EnergyResult SubEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd, int st,
      int en, bool build_structure = false) const;
  EnergyResult TotalEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd,
      bool build_structure = false) const;

  bool IsValid(std::string* reason = nullptr) const {
    CHECK_COND(multiloop_c == ZERO_E, "multiloop_c must be zero");
    return base::ModelIsValid(*this, reason);
  }

  void LoadPseudofreeEnergy(std::vector<Energy> pf_paired, std::vector<Energy> pf_unpaired) {
    pf.Load(std::move(pf_paired), std::move(pf_unpaired));
  }

  void LoadFromModelPath(const std::string& path);
  void LoadRandom(std::mt19937& eng);

 private:
  friend class ModelMixin<Model>;

  // This is private to prevent construction on the stack, since this structure is large.
  Model() = default;

  Energy SubEnergyInternal(const Primary& r, const Secondary& s, int st, int en, int stack_st,
      int stack_en, bool use_given_ctds, Ctds* ctd, std::unique_ptr<Structure>* struc) const;
};

}  // namespace mrna::md::stack

#endif  // BACKENDS_STACK_ENERGY_MODEL_H_
