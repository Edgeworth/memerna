// Copyright 2021 E.
#include "compute/subopt/brute.h"

#include <set>
#include <stack>
#include <utility>
#include <vector>

#include "compute/brute.h"
#include "compute/energy/structure.h"
#include "util/splaymap.h"

namespace mrna::subopt {

std::vector<subopt::SuboptResult> SuboptimalBruteForce(
    const Primary& r, const energy::EnergyModel& em, int max_structures) {
  auto res = BruteForce().Run(r, em, max_structures, false, false);
  return std::vector<subopt::SuboptResult>(res.subopts.begin(), res.subopts.end());
}

}  // namespace mrna::subopt
