#include "fold/fold.h"

namespace memerna {
namespace fold {

namespace {

energy_t best;
std::vector<int> best_p;
std::vector<std::pair<int, int>> base_pairs;

void FoldBruteForceInternal(int idx) {
  if (idx == int(base_pairs.size())) {
    energy_t e = energy::ComputeEnergy(nullptr);
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

folded_rna_t FoldBruteForce(const rna_t& rna) {
  SetRna(rna);
  p = std::vector<int>(r.size(), -1);
  best = constants::MAX_E;
  base_pairs.clear();
  for (int st = 0; st < int(r.size()); ++st) {
    for (int en = st + constants::HAIRPIN_MIN_SZ + 1; en < int(r.size()); ++en) {
      if (CanPair(r[st], r[en]) && IsNotLonely(st, en))
        base_pairs.emplace_back(st, en);
    }
  }
  FoldBruteForceInternal(0);
  p = best_p;
  return {rna, p, energy::ComputeEnergy()};
}

}
}
