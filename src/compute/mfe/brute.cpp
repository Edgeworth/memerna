// Copyright 2016 E.
#include "compute/mfe/brute.h"

#include <set>
#include <stack>
#include <utility>
#include <vector>

#include "compute/energy/structure.h"
#include "compute/mfe/mfe.h"
#include "compute/subopt/brute.h"
#include "util/splaymap.h"

namespace mrna::mfe {

subopt::SuboptResult MfeBruteForce(Primary r, const energy::EnergyModel& em) {
  return std::move(subopt::SuboptimalBruteForce(std::move(r), em, 1)[0]);
}

}  // namespace mrna::mfe
