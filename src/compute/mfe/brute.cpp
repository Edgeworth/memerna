// Copyright 2016 Eliot Courtney.
#include "compute/mfe/brute.h"

#include <utility>
#include <vector>

#include "compute/subopt/brute.h"
#include "model/primary.h"

namespace mrna::mfe {

subopt::SuboptResult MfeBruteForce(Primary r, const energy::EnergyModel& em) {
  return std::move(subopt::SuboptimalBruteForce(std::move(r), em, 1)[0]);
}

}  // namespace mrna::mfe
