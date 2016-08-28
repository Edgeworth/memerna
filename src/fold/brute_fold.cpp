#include "fold/fold.h"
#include "parsing.h"

namespace memerna {
namespace fold {

namespace {

energy_t best;
std::vector<int> best_p;
std::vector<std::pair<int, int>> base_pairs;

void FoldBruteForceInternal(int idx) {
  if (idx == int(base_pairs.size())) {
    energy_t e = energy::ComputeEnergy();
    if (e < best) {
      best = e;
      best_p = p;
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
    if (p[i] != -1) {
      can_take = false;
      break;
    }
  }
  if (can_take) {
    p[pair.first] = pair.second;
    p[pair.second] = pair.first;
    FoldBruteForceInternal(idx + 1);
    p[pair.first] = -1;
    p[pair.second] = -1;
  }
}

}

folded_rna_t FoldBruteForce(const rna_t& rna, fold_state_t*) {
  SetRna(rna);
  p = std::vector<int>(r.size(), -1);
  best = constants::MAX_E;
  best_p.clear();
  base_pairs.clear();
  // Add base pairs in order of increasing st, then en.
  for (int st = 0; st < int(r.size()); ++st) {
    for (int en = st + constants::HAIRPIN_MIN_SZ + 1; en < int(r.size()); ++en) {
      if (CanPair(r[st], r[en]) && IsNotLonely(st, en))
        base_pairs.emplace_back(st, en);
    }
  }
  FoldBruteForceInternal(0);
  p = best_p;
  folded_rna_t frna = {rna, best_p, energy::ComputeEnergy()};
  return frna;
}

}
}
