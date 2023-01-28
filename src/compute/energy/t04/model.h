// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_ENERGY_T04_MODEL_H_
#define COMPUTE_ENERGY_T04_MODEL_H_

#include <memory>

#include "compute/energy/common/model.h"
#include "compute/energy/common/t04like/model_mixin.h"
#include "compute/energy/energy.h"
#include "compute/energy/structure.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"

namespace mrna::erg::t04 {

class Model : public ModelMixin<Model>, public T04ModelMixin {
 public:
  EnergyResult SubEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd, int st,
      int en, bool build_structure = false) const;
  EnergyResult TotalEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd,
      bool build_structure = false) const;

 private:
  friend class ModelMixin<Model>;

  // This is private to prevent construction on the stack, since this structure is large.
  Model() = default;

  Energy SubEnergyInternal(const Primary& r, const Secondary& s, int st, int en,
      bool use_given_ctds, Ctds* ctd, std::unique_ptr<Structure>* struc) const;
};

}  // namespace mrna::erg::t04

#endif  // COMPUTE_ENERGY_T04_MODEL_H_
