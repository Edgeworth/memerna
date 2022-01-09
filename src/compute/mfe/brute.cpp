// Copyright 2016 Eliot Courtney.
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

Computed MfeBruteForce(const Primary& r, const energy::EnergyModel& em) {
  return subopt::SuboptimalBruteForce(r, em, 1)[0];
}

}  // namespace mrna::mfe
