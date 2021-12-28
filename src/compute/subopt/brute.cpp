// Copyright 2021 E.
#include "compute/subopt/brute.h"

#include <set>
#include <stack>
#include <utility>
#include <vector>

#include "compute/brute.h"
#include "compute/energy/globals.h"
#include "compute/energy/structure.h"
#include "compute/mfe/globals.h"
#include "model/parsing.h"
#include "util/splaymap.h"

namespace mrna::subopt {

std::vector<Computed> SuboptimalBruteForce(
    const Primary& r, const energy::EnergyModel& em, int max_structures_) {
  auto res = BruteForce().Run(r, em, max_structures_, false, false);
  return std::vector<Computed>(res.best_computeds.begin(), res.best_computeds.end());
}

}  // namespace mrna::subopt
