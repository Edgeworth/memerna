// Copyright 2021 E.
#include "compute/subopt/brute.h"

#include <set>
#include <utility>
#include <vector>

#include "compute/brute.h"
#include "compute/subopt/subopt.h"
#include "model/primary.h"

namespace mrna::subopt {

std::vector<subopt::SuboptResult> SuboptimalBruteForce(
    Primary r, const energy::EnergyModel& em, int max_structures) {
  auto res = BruteForce().Run(std::move(r), em, max_structures, false, false);
  return std::vector<subopt::SuboptResult>(res.subopts.begin(), res.subopts.end());
}

}  // namespace mrna::subopt
