// Copyright 2022 Eliot Courtney.
#ifndef COMPUTE_ENERGY_T22_MODEL_H_
#define COMPUTE_ENERGY_T22_MODEL_H_

#include <fmt/core.h>

#include <cassert>
#include <memory>
#include <random>
#include <string>

#include "compute/energy/common/model.h"
#include "compute/energy/common/t04like/model_mixin.h"
#include "compute/energy/energy.h"
#include "compute/energy/structure.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"

namespace mrna::erg::t22 {

// TODO(0): Implement.
class Model : public ModelMixin<Model>, public T04ModelMixin {
 public:
  Energy penultimate_stack[4][4][4][4] = {};

  EnergyResult SubEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd, int st,
      int en, bool build_structure = false) const;
  EnergyResult TotalEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd,
      bool build_structure = false) const;

  // Computes the penultimate penalty for the closing base pair at (st, en),
  // assuming that (st-1, en+1) is paired. This computes the penultimate penalty
  // facing inwards. To compute for facing outwards, reverse st and en.
  [[nodiscard]] constexpr Energy PenultimatePenalty(const Primary& r, int st, int en) const {
    assert(st > 0 && en < static_cast<int>(r.size()) - 1);
    return penultimate_stack[r[st - 1]][r[st]][r[en]][r[en + 1]];
  }

  // Computes the penalty for a stack of the given length, ending at (ist, ien).
  [[nodiscard]] constexpr Energy StackPenalty(
      const Primary& r, int ist, int ien, int len, std::unique_ptr<Structure>* s = nullptr) const {
    assert(len >= 0);
    if (len <= 1) return ZERO_E;
    int ost = ist - len + 1;
    int oen = ien + len - 1;
    auto inner = PenultimatePenalty(r, ist, ien);
    auto outer = PenultimatePenalty(r, oen, ost);
    if (s) {
      if (s) (*s)->AddNote("{}e - inner penultimate penalty at ({}, {})", inner, ist, ien);
      if (s) (*s)->AddNote("{}e - outer penultimate penalty at ({}, {})", outer, ost, oen);
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

  Energy SubEnergyInternal(const Primary& r, const Secondary& s, int st, int en, int stack_len,
      bool use_given_ctds, Ctds* ctd, std::unique_ptr<Structure>* struc) const;
};

}  // namespace mrna::erg::t22

#endif  // COMPUTE_ENERGY_T22_MODEL_H_
