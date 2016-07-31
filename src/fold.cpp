#include <queue>
#include "fold.h"
#include "globals.h"

namespace memerna {
namespace fold {

namespace {

#define UPDATE_CACHE(arr, us_st, us_en, value) \
  do { \
    energy_t macro_upd_value_ = (value); \
    if (macro_upd_value_ < CAP_E && macro_upd_value_ < arr[us_st][us_en]) { \
      /*printf("Updating " #arr " %d %d %d => %d\n  " #value "\n", us_st, us_en, arr[us_st][us_en], value);*/ \
      arr[us_st][us_en] = macro_upd_value_; \
    } \
  } while (0)

energy_t FoldInternal() {
  int N = int(r.size());
  assert(N > 0);
  // Automatically initialised to MAX_E.
  array2d_t<energy_t> paired(r.size() + 1), unpaired(r.size() + 1);
  assert(paired[N - 1][0] == MAX_E);

  // TODO: check bounds
  // TODO: au/gu penalty?
  // sz includes st and en.
  for (int sz = HAIRPIN_MIN_SZ + 2; sz <= N; ++sz) {
    int en = sz - 1;
    for (int st = 0; st < N - sz + 1; ++st) {
      // Update paired.
      // Internal loops, bulge loops, and stacking. TODO: Lyngso's
      // TODO: can reduce the range of these loops.
      if (CanPair(r[st], r[en])) {
        for (int ist = st + 1; ist < en - 1; ++ist) {
          for (int ien = ist + 1; ien < en; ++ien) {
            UPDATE_CACHE(paired, sz, st, energy::TwoLoop(st, en, ist, ien) + paired[ien - ist + 1][ist]);
          }
        }
      }
      // Hairpin loops.
      UPDATE_CACHE(paired, sz, st, energy::HairpinEnergy(st, en));

      // Multiloops. Look at range [st + 1, en - 1].
      // TODO: can reduce range
      energy_t augu_penalty = IsAuGu(r[st], r[en]) ? energy::AUGU_PENALTY * 0 : 0;
      for (int pivsz = HAIRPIN_MIN_SZ + 2; pivsz < sz - 1; ++pivsz) {
        // Include AU/GU penalty for ending multiloop helix.
        UPDATE_CACHE(paired, sz, st, unpaired[pivsz][st + 1] +
            unpaired[sz - pivsz - 2][st + 1 + pivsz] + augu_penalty + multiloop_hack_a);
      }

      // Update unpaired.
      // Choose |st| to be unpaired.
      UPDATE_CACHE(unpaired, sz, st, unpaired[sz - 1][st + 1] + multiloop_hack_b);
      // Pair here.
      UPDATE_CACHE(unpaired, sz, st, paired[sz][st]);
      // TODO: can reduce range, probably
      for (int pivsz = HAIRPIN_MIN_SZ + 2; pivsz < sz; ++pivsz) {
        int rem = sz - pivsz;
        augu_penalty = IsAuGu(r[st], r[st + pivsz - 1]) ? energy::AUGU_PENALTY * 0 : 0;
        // Split.
        UPDATE_CACHE(unpaired, sz, st, paired[pivsz][st] + unpaired[rem][st + pivsz] + augu_penalty);
        // Choose rest to be empty.
        UPDATE_CACHE(unpaired, sz, st, paired[pivsz][st] + multiloop_hack_b * rem + augu_penalty);
      }

      en++;
    }
  }

  // Exterior loop calculation. There can be no paired base on exterior[en].
  std::vector<energy_t> exterior(std::size_t(N + 1), 0);
  std::vector<int> exterior_sts(std::size_t(N + 1), -1);
  for (int en = 1; en <= N; ++en) {
    exterior[en] = exterior[en - 1];
    for (int st = 0; st < en; ++st) {
      energy_t augu_penalty = IsAuGu(r[st], r[en - 1]) ? energy::AUGU_PENALTY * 0 : 0;
      energy_t value = exterior[st] + paired[en - st][st] + augu_penalty;
      if (value < exterior[en]) {
        exterior[en] = value;
        exterior_sts[en] = st;
      }
    }
  }

  // Backtrace.
  p = std::vector<int>(r.size(), -1);
  // sz, st, paired
  std::queue<std::tuple<int, int, bool>> q;
  int ext_st = N;
  while (ext_st > 0) {
    int nst = exterior_sts[ext_st];
    printf("%d %d\n", ext_st, nst);
    if (nst == -1) {
      --ext_st;
      continue;
    }
    q.emplace(ext_st - nst, nst, true);
    p[nst] = ext_st - 1;
    p[ext_st - 1] = nst;
    ext_st = nst;
  }

  while (!q.empty()) {
    int sz, st;
    bool is_paired;
    std::tie(sz, st, is_paired) = q.front();
    int en = st + sz - 1;
    q.pop();
    if (is_paired) {
      // TODO
    } else {
      // TODO
    }
  }
  return exterior[N];
}

#undef UPDATE_CACHE

}

energy_t Fold(const rna_t& rna, std::unique_ptr<structure::Structure>* s) {
  r = rna;
  return FoldInternal();
}

namespace {

energy_t best;
std::vector<int> best_p;
std::vector<std::pair<int, int>> base_pairs;

void FoldBruteForceInternal(int idx) {
  // Reached the end. r and p already set, so directly call ComputeEnergy
  if (idx == int(base_pairs.size())) {
    energy_t e = energy::ComputeEnergy(0, int(r.size() - 1), nullptr);
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

energy_t FoldBruteForce(const rna_t& rna, std::unique_ptr<structure::Structure>* s) {
  r = rna;
  p = std::vector<int>(r.size(), -1);
  best = MAX_E;
  base_pairs.clear();
  for (int st = 0; st < int(r.size()); ++st) {
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < int(r.size()); ++en) {
      if (CanPair(r[st], r[en]))
        base_pairs.emplace_back(st, en);
    }
  }
  FoldBruteForceInternal(0);
  p = best_p;
  return energy::ComputeEnergy(0, int(r.size() - 1), s);
}

}
}
