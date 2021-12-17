// Copyright 2016 E.
#include "compute/mfe/brute.h"

#include <set>
#include <stack>
#include <utility>
#include <vector>

#include "compute/energy/globals.h"
#include "compute/energy/structure.h"
#include "compute/mfe/globals.h"
#include "compute/mfe/mfe.h"
#include "compute/subopt/brute.h"
#include "model/parsing.h"
#include "util/splaymap.h"

namespace mrna::mfe {

computed_t MfeBruteForce(const primary_t& r, const energy::EnergyModel& em) {
  return subopt::SuboptimalBruteForce(r, em, 1)[0];
}

}  // namespace mrna::mfe
