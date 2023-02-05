#ifndef API_BRUTE_H_
#define API_BRUTE_H_

namespace mrna {

struct BruteResult {
  // Suboptimal structures:
  // TODO(2): use splayset here?
  std::multiset<subopt::SuboptResult, SuboptCmp> subopts;

  // Partition function:
  part::Part part;
  BoltzProbs prob;
};

}

#endif  // API_BRUTE_H_
