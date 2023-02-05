// Copyright 2022 Eliot Courtney.
#ifndef MODELS_T22_ENERGY_MODEL_H_
#define MODELS_T22_ENERGY_MODEL_H_

#include <fmt/core.h>

#include <cassert>
#include <memory>
#include <random>
#include <string>

#include "api/energy/energy.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "model/structure.h"
#include "models/common/model.h"
#include "models/t04/energy/model_mixin.h"

namespace mrna::md::t22::erg {

using mrna::erg::EnergyResult;

class Model : public ModelMixin<Model>, public t04::erg::T04ModelMixin {
 public:
  Energy penultimate_stack[4][4][4][4] = {};

  EnergyResult SubEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd, int st,
      int en, bool build_structure = false) const;
  EnergyResult TotalEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd,
      bool build_structure = false) const;

  // Computes the penalty for a stack of the given length, ending at (ist, ien).
  // Handles bulge loops.
  [[nodiscard]] constexpr Energy StackPenalty(const Primary& r, const Secondary& s, int ost,
      int oen, int ist, int ien, std::unique_ptr<Structure>* struc = nullptr) const {
    assert(ost >= 0 && oen < static_cast<int>(r.size()));
    assert(ist > 0 && ien < static_cast<int>(r.size()) - 1);
    assert(ist > ost && ien < oen);
    // Check for single unpaired bases to treat as continuous.
    int ost_next = s[ost + 1] == -1 ? ost + 2 : ost + 1;
    int oen_prev = s[oen - 1] == -1 ? oen - 2 : oen - 1;
    assert(ost_next - ost + oen - oen_prev < 4);

    int ist_prev = s[ist - 1] == -1 ? ist - 2 : ist - 1;
    int ien_next = s[ien + 1] == -1 ? ien + 2 : ien + 1;
    assert(ist - ist_prev + ien_next - ien < 4);

    auto inner = penultimate_stack[r[ist_prev]][r[ist]][r[ien]][r[ien_next]];
    auto outer = penultimate_stack[r[oen_prev]][r[oen]][r[ost]][r[ost_next]];
    if (struc) {
      (*struc)->AddNote("{}e - inner penultimate penalty at ({}, {})", inner, ist, ien);
      (*struc)->AddNote("{}e - outer penultimate penalty at ({}, {})", outer, ost, oen);
    }

    return inner + outer;
  }

 protected:
  void LoadFromDir(const std::string& data_dir);

  void LoadRandom(std::mt19937& eng);

 private:
  friend class ModelMixin<Model>;

  // This is private to prevent construction on the stack, since this structure is large.
  Model() = default;

  Energy SubEnergyInternal(const Primary& r, const Secondary& s, int st, int en, int stack_st,
      int stack_en, bool use_given_ctds, Ctds* ctd, std::unique_ptr<Structure>* struc) const;
};

}  // namespace mrna::md::t22::erg

#endif  // MODELS_T22_ENERGY_MODEL_H_
