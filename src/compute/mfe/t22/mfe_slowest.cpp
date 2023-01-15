// Copyright 2023 Eliot Courtney.
#include <algorithm>
#include <memory>

#include "compute/dp.h"
#include "compute/energy/t22/model.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/energy.h"
#include "model/primary.h"
#include "util/array.h"

namespace mrna::mfe::t22 {

// TODO(0): Implement.
DpArray MfeSlowest(const Primary& r, const erg::t22::ModelPtr& em) {
  bug();
  return DpArray::empty();
}

}  // namespace mrna::mfe::t22
