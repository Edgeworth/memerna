// Copyright 2016 Eliot Courtney.
#include "model/branch.h"

#include <algorithm>
#include <cassert>
#include <vector>

#include "model/base.h"

namespace mrna {

std::vector<int> GetBranchCounts(const Secondary& s) {
  std::vector<int> branch_count(s.size(), 0);
  std::vector<int> q;
  for (int i = 0; i < static_cast<int>(s.size()); ++i) {
    if (s[i] == -1) continue;
    if (s[i] > i) {
      // Exterior loop counts a multiloop for CTDs.
      if (q.empty()) branch_count[i] = 2;
      q.push_back(i);

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
      q.pop_back();
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

int MaxNumContiguous(const Primary& r) {
  int num_contig = 0;
  int max_num_contig = 0;
  Base prev = -1;
  for (auto b : r) {
    if (b == prev)
      num_contig++;
    else
      num_contig = 1;
    prev = b;
    max_num_contig = std::max(max_num_contig, num_contig);
  }
  return max_num_contig;
}

}  // namespace mrna
