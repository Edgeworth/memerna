// Copyright 2016 Eliot Courtney.
#ifndef MODELS_T04_ENERGY_MODEL_H_
#define MODELS_T04_ENERGY_MODEL_H_

#include <memory>

#include "api/energy/energy.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "model/structure.h"
#include "models/common/model.h"
#include "models/t04/energy/model_mixin.h"

namespace mrna::md::t04 {

using mrna::erg::EnergyResult;

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

}  // namespace mrna::md::t04

#endif  // MODELS_T04_ENERGY_MODEL_H_
