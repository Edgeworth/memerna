// Copyright 2023 Eliot Courtney.
#include "compute/dp.h"
#include "compute/energy/t22/model.h"
#include "model/energy.h"
#include "model/primary.h"
#include "util/error.h"

namespace mrna::mfe::t22 {

// TODO(0): Implement. Any way to generalise/pull out code?
ExtArray MfeExterior(const Primary& r, const erg::t22::Model::Ptr& /*em*/, const DpArray& /*dp*/) {
  auto ext = ExtArray(r.size() + 1, MAX_E);
  bug();
  return ext;
}

}  // namespace mrna::mfe::t22
