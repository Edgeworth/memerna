// Copyright 2022 Eliot Courtney.
#ifndef COMPUTE_ENERGY_T22_MODEL_H_
#define COMPUTE_ENERGY_T22_MODEL_H_

#include <cassert>
#include <cmath>
#include <cstdint>
#include <deque>
#include <memory>
#include <string>
#include <unordered_map>

#include "compute/energy/common/model.h"
#include "compute/energy/common/t04like/model_mixin.h"
#include "compute/energy/energy.h"
#include "compute/energy/energy_cfg.h"
#include "compute/energy/structure.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/argparse.h"

namespace mrna::erg::t22 {

// TODO(0): Implement.
class Model : public ModelMixin<Model>, public T04ModelMixin {
 public:
  Energy penultimate_stack[4][4][4][4] = {};

  Energy Hairpin(const Primary& r, int st, int en, bool prev_paired,
      std::unique_ptr<Structure>* s = nullptr) const;
  Energy Bulge(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const;
  Energy InternalLoop(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const;
  Energy TwoLoop(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const;
  Energy MultiloopEnergy(const Primary& r, const Secondary& s, int st, int en,
      std::deque<int>* branches, bool use_given_ctds, Ctds* ctd,
      std::unique_ptr<Structure>* sstruc = nullptr) const;

  EnergyResult SubEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd, int st,
      int en, bool build_structure = false) const;
  EnergyResult TotalEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd,
      bool build_structure = false) const;

  // Computes the penultimate penalty for the closing base pair at (st, en),
  // assuming that (st-1, en+1) is paired. This computes the penultimate penalty
  // facing inwards. To compute for facing outwards, reverse st and en.
  Energy PenultimatePenalty(const Primary& r, int st, int en) const {
    assert(st > 0 && en < static_cast<int>(r.size()) - 1);
    return penultimate_stack[r[st - 1]][r[st]][r[en]][r[en + 1]];
  }

  // Computes the penultimate penalty for the closing base pair at (st, en),
  // if (st-1, en+1) is paired.
  Energy PenultimatePenalty(const Primary& r, const Secondary& s, int st, int en) const {
    // No inner penultimate penalty if there can be no previous base pair..
    if (st == 0 || en == static_cast<int>(r.size()) - 1) return ZERO_E;
    // No inner penultimate penalty if there is no previous base pair.
    if (s[st - 1] != en + 1) return ZERO_E;
    return PenultimatePenalty(r, st, en);
  }

 protected:
  void LoadFromDir(const std::string& data_dir);

  void LoadRandom(std::mt19937& eng);

 private:
  friend class ModelMixin<Model>;

  // This is private to prevent construction on the stack, since this structure is large.
  Model() = default;
};

}  // namespace mrna::erg::t22

#endif  // COMPUTE_ENERGY_T22_MODEL_H_
