// Copyright 2023 Eliot Courtney.
#include "compute/traceback/t22/traceback.h"

#include "compute/dp.h"
#include "model/primary.h"
#include "util/error.h"

namespace mrna::tb::t22 {

// TODO(0): Implement. Think if can generalise this.
TracebackResult Traceback(const Primary& /*r*/, const erg::t22::Model::Ptr& /*em*/,
    const DpArray& /*dp*/, const ExtArray& /*ext*/) {
  bug();
  return TracebackResult{};
}

}  // namespace mrna::tb::t22
