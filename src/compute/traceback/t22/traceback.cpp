// Copyright 2023 Eliot Courtney.
#include "compute/traceback/t22/traceback.h"

#include <algorithm>
#include <memory>
#include <stack>

#include "compute/dp.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/array.h"

namespace mrna::tb::t22 {

// TODO(0): Implement. Think if can generalise this.
TracebackResult Traceback(
    const Primary& r, const erg::t22::ModelPtr& em, const DpArray& dp, const ExtArray& ext) {
  bug();
  return TracebackResult{};
}

}  // namespace mrna::tb::t22
