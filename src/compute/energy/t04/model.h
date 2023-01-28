// Copyright 2016 E.
#ifndef COMPUTE_ENERGY_T04_MODEL_H_
#define COMPUTE_ENERGY_T04_MODEL_H_

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

namespace mrna::erg::t04 {

class Model : public ModelMixin<Model>, public T04ModelMixin {
 public:
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

  EnergyResult SubEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd, int st,
      int en, bool build_structure = false) const;
  EnergyResult TotalEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd,
      bool build_structure = false) const;

 private:
  friend class ModelMixin<Model>;

  // This is private to prevent construction on the stack, since this structure is large.
  Model() = default;
};

}  // namespace mrna::erg::t04

#endif  // COMPUTE_ENERGY_T04_MODEL_H_
