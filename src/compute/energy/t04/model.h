// Copyright 2016 Eliot Courtney.
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
 private:
  friend class ModelMixin<Model>;

  // This is private to prevent construction on the stack, since this structure is large.
  Model() = default;
};

}  // namespace mrna::erg::t04

#endif  // COMPUTE_ENERGY_T04_MODEL_H_
