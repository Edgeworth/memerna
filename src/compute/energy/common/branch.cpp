// Copyright 2016 Eliot Courtney.
#include "compute/energy/common/branch.h"

#include <cassert>
#include <memory>
#include <stack>
#include <vector>

namespace mrna::erg {

std::vector<int> GetBranchCounts(const Secondary& s) {
  std::vector<int> branch_count(s.size(), 0);
  std::stack<int> q;
  for (int i = 0; i < static_cast<int>(s.size()); ++i) {
    if (s[i] == -1) continue;
    if (s[i] > i) {
      // Exterior loop counts a multiloop for CTDs.
      if (q.empty()) branch_count[i] = 2;
      q.push(i);

      // Look at all the children.
      int count = 0;
      for (int j = i + 1; j < s[i]; ++j) {
        if (s[j] != -1) {
          j = s[j];
          count++;
        }
      }
      for (int j = i + 1; j < s[i]; ++j) {
        if (s[j] != -1) {
          branch_count[j] = count;
          j = s[j];
        }
      }
      branch_count[s[i]] = count;
    } else {
      q.pop();
    }
  }
  return branch_count;
}

void AddBranchCtdsToBaseCtds(
    const std::deque<int>& branches, const BranchCtd& branch_ctd, Ctds* ctd) {
  assert(branches.size() == branch_ctd.size());
  for (int i = 0; i < static_cast<int>(branches.size()); ++i) {
    // Only write it into one side. If it's for an outer loop, it will be the right side, since we
    // swap the indices in that case.
    (*ctd)[branches[i]] = branch_ctd[i].first;
  }
}

}  // namespace mrna::erg
