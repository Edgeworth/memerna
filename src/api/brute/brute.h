#ifndef API_BRUTE_H_
#define API_BRUTE_H_

#include <set>

#include "api/part.h"
#include "api/subopt/subopt.h"

namespace mrna::brute {

struct SuboptCmp {
  bool operator()(const subopt::SuboptResult& a, const subopt::SuboptResult& b) const {
    // Kept in a multiset, so this is just used for ordering, not deduplication.
    // There should be no duplicates added anyway. Single comparison to keep it fast.
    return a.energy < b.energy;
  }
};

struct BruteResult {
  // Suboptimal structures:
  // TODO(2): use splayset here?
  std::multiset<subopt::SuboptResult, SuboptCmp> subopts;

  // Partition function:
  part::Part part;
  part::BoltzProbs prob;
};

}  // namespace mrna::brute

#endif  // API_BRUTE_H_
