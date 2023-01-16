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

// TODO(0): Implement. Any way to generalise/pull out code?
DpArray MfeSlowest(const Primary& r, const erg::t22::Model::Ptr& em) {
  auto dp = DpArray(r.size() + 1, MAX_E);
  bug();
  return dp;
}

}  // namespace mrna::mfe::t22
