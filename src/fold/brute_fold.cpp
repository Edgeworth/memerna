#include <set>
#include "fold/brute_fold.h"
#include "fold/globals.h"
#include "fold/fold_internal.h"
#include "parsing.h"
#include "energy/structure.h"

namespace memerna {
namespace fold {

namespace internal {

std::vector<int> GetBranchCounts(const std::vector<int>& p) {
  std::vector<int> branch_count(p.size(), 0);
  std::stack<int> q;
  for (int i = 0; i < int(p.size()); ++i) {
    if (p[i] == -1) continue;
    if (p[i] > i) {
      // Exterior loop counts a multiloop for CTDs.
      if (q.empty())
        branch_count[i] = 2;
      q.push(i);

      // Look at all the children.
      int count = 0;
      for (int j = i + 1; j < p[i]; ++j) {
        if (p[j] != -1) {
          j = p[j];
          count++;
        }
      }
      for (int j = i + 1; j < p[i]; ++j) {
        if (p[j] != -1) {
          branch_count[j] = count;
          j = p[j];
        }
      }
    } else {
      q.pop();
    }
  }
  return branch_count;
}

}

namespace {

using namespace internal;

std::set<computed_t> best_computeds;
std::vector<std::pair<int, int>> base_pairs;
std::vector<int> branch_count;
int max_structures;

void AddAllCombinations(int idx) {
  // Base case
  if (idx == int(gr.size())) {
    printf("%s\n%s\n", parsing::PrimaryToString(gr).c_str(), parsing::PairsToDotBracket(gp).c_str());
    for (auto ctd : gctd) {
      printf("%s, ", energy::CtdToName(ctd));
    }
    printf("\n");
    for (auto cnt : branch_count)
      printf("%d, ", cnt);
    printf("\n");
    const auto computed = energy::ComputeEnergyWithCtds({{gr, gp}, gctd, 0}, gem);
    if (int(best_computeds.size()) < max_structures || best_computeds.rbegin()->energy > computed.energy)
      best_computeds.insert(std::move(computed));
    if (int(best_computeds.size()) > max_structures)
      best_computeds.erase(--best_computeds.end());
    return;
  }
  // If we already set this, this isn't a valid base pair, it's not part of a multiloop, or it's the right side
  // of a branch, can't or don't want to set ctds so continue.
  if (gctd[idx] != CTD_NA || gp[idx] == -1 || branch_count[idx] < 2 || gp[idx] < idx) {
    AddAllCombinations(idx + 1);
    return;
  }
  const int N = int(gr.size());
  const bool lu_exists = idx - 1 >= 0 && gp[idx - 1] == -1;
  const bool lu_shared = lu_exists && idx - 2 >= 0 && gp[idx - 2] != -1;
  const bool lu_usable = lu_exists && (!lu_shared || (
      gctd[gp[idx - 2]] != CTD_3_DANGLE &&
          gctd[gp[idx - 2]] != CTD_TERMINAL_MISMATCH &&
          gctd[gp[idx - 2]] != CTD_RIGHT_MISMATCH_COAX_WITH_PREV));
  const bool ru_exists = gp[idx] + 1 < N && gp[gp[idx] + 1] == -1;
  const bool ru_shared = ru_exists && gp[idx] + 2 < N && gp[gp[idx] + 2] != -1;
  const bool ru_usable = ru_exists && (!ru_shared || (
      gctd[gp[idx] + 2] != CTD_5_DANGLE &&
          gctd[gp[idx] + 2] != CTD_TERMINAL_MISMATCH &&
          gctd[gp[idx] + 2] != CTD_LEFT_MISMATCH_COAX_WITH_NEXT));
  // Even if the next branch is an outer branch, everything will be magically handled.
  // CTD_UNUSED
  gctd[idx] = CTD_UNUSED;
  AddAllCombinations(idx + 1);

  // CTD_3_DANGLE
  if (ru_usable) {
    gctd[idx] = CTD_3_DANGLE;
    AddAllCombinations(idx + 1);
  }

  // CTD_5_DANGLE
  if (lu_usable) {
    gctd[idx] = CTD_5_DANGLE;
    AddAllCombinations(idx + 1);
  }

  // CTD_TERMINAL_MISMATCH
  if (ru_usable && lu_usable) {
    gctd[idx] = CTD_TERMINAL_MISMATCH;
    AddAllCombinations(idx + 1);
  }

  // Check that the next branch hasn't been set already.
  // CTD_LEFT_MISMATCH_COAX_WITH_NEXT
  if (lu_usable && ru_usable && ru_shared && gctd[gp[idx] + 2] == CTD_NA) {
    gctd[idx] = CTD_LEFT_MISMATCH_COAX_WITH_NEXT;
    gctd[gp[idx] + 2] = CTD_LEFT_MISMATCH_COAX_WITH_PREV;
    AddAllCombinations(idx + 1);
  }

  // Check that the previous branch hasn't been set already.
  // CTD_RIGHT_MISMATCH_COAX_WITH_NEXT
  if (lu_usable && lu_shared && ru_usable && gctd[gp[idx - 2]] == CTD_NA) {
    gctd[gp[idx - 2]] = CTD_RIGHT_MISMATCH_COAX_WITH_NEXT;
    gctd[idx] = CTD_RIGHT_MISMATCH_COAX_WITH_PREV;
    AddAllCombinations(idx + 1);
  }

  // CTD_FLUSH_COAX_WITH_NEXT
  if (gp[idx] + 1 < N && gp[gp[idx] + 1] != -1 && gctd[gp[idx] + 1] == CTD_NA) {
    gctd[idx] = CTD_FLUSH_COAX_WITH_NEXT;
    gctd[gp[idx] + 1] = CTD_FLUSH_COAX_WITH_PREV;
    AddAllCombinations(idx + 1);
  }

  // Reset back to NA.
  gctd[idx] = CTD_NA;
}

void FoldBruteForceInternal(int idx) {
  if (idx == int(base_pairs.size())) {
    // Small optimisation for case when we're just getting one structure.
    if (max_structures == 1) {
      auto computed = energy::ComputeEnergy({gr, gp}, gem);
      if (best_computeds.empty() || computed.energy < best_computeds.begin()->energy)
        best_computeds.insert(std::move(computed));
      if (best_computeds.size() == 2)
        best_computeds.erase(--best_computeds.end());
    } else {
      // Precompute whether things are multiloops or not.
      branch_count = internal::GetBranchCounts(gp);
      AddAllCombinations(0);
    }

    return;
  }
  // Don't take this base pair.
  FoldBruteForceInternal(idx + 1);

  // Take this base pair.
  bool can_take = true;
  const auto& pair = base_pairs[idx];
  // Only need to check in the range of this base pair. Since we ordered by
  // increasing st, anything at or after this will either be the start of something starting at st, or
  // something ending, both of which conflict with this base pair.
  for (int i = pair.first; i <= pair.second; ++i) {
    if (gp[i] != -1) {
      can_take = false;
      break;
    }
  }
  if (can_take) {
    gp[pair.first] = pair.second;
    gp[pair.second] = pair.first;
    FoldBruteForceInternal(idx + 1);
    gp[pair.first] = -1;
    gp[pair.second] = -1;
  }
}

}

std::vector<computed_t> FoldBruteForce(const primary_t& r,
    const energy::EnergyModel& em, int max_structures_) {
  internal::SetGlobalState(r, em);
  best_computeds.clear();
  base_pairs.clear();
  max_structures = max_structures_;
  // Add base pairs in order of increasing st, then en.
  for (int st = 0; st < int(r.size()); ++st) {
    for (int en = st + constants::HAIRPIN_MIN_SZ + 1; en < int(r.size()); ++en) {
      if (internal::ViableFoldingPair(st, en))
        base_pairs.emplace_back(st, en);
    }
  }
  FoldBruteForceInternal(0);
  return std::vector<computed_t>(best_computeds.begin(), best_computeds.end());
}

}
}
