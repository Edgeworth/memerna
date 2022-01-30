// Copyright 2021 E.
#include "compute/subopt/brute.h"

#include <set>
#include <utility>
#include <vector>

#include "compute/brute/brute.h"
#include "compute/subopt/subopt.h"
#include "model/primary.h"

namespace mrna::subopt {

std::vector<subopt::SuboptResult> SuboptimalBruteForce(
    Primary r, const energy::EnergyModel& em, int strucs) {
  auto res = brute::BruteForce().Run(std::move(r), em, strucs, false, false);
  return std::vector<subopt::SuboptResult>(res.subopts.begin(), res.subopts.end());
}

}  // namespace mrna::subopt
