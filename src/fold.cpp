#include <queue>
#include "fold.h"

namespace memerna {
namespace fold {

#define UPDATE_CACHE(a, sz, st, value) \
  do { \
    energy_t macro_upd_value_ = (value); \
    if (macro_upd_value_ < CAP_E && macro_upd_value_ < arr[sz][st][a]) { \
      /*printf("Updating " #arr " %d %d %d => %d\n  " #value "\n", sz, st, arr[sz][st][a], value);*/ \
      arr[sz][st][a] = macro_upd_value_; \
    } \
  } while (0)

enum {
  P,
  U,
  U_ST,
  U_USED,
  U_RCOAX,
  ARR_SIZE
};

energy_t Fold(std::unique_ptr<structure::Structure>* s) {
  int N = int(r.size());
  assert(N > 0);
  // Automatically initialised to MAX_E.
  array3d_t<energy_t, ARR_SIZE> arr(r.size() + 1);
  assert(arr[N - 1][0][U_ST] == MAX_E);

  // TODO: check bounds
  // TODO: au/gu penalty?
  // sz includes st and en.
  static_assert(HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");
  for (int sz = HAIRPIN_MIN_SZ + 2; sz <= N; ++sz) {
    int en = sz - 1;
    for (int st = 0; st < N - sz + 1; ++st) {
      base_t stb = r[st], st1b = r[st + 1], st2b = r[st + 2], enb = r[en];

      // Update paired - only if can actually pair.
      if (CanPair(r[st], r[en])) {
        // Internal loops, bulge loops, and stacking. TODO: Lyngso's
        // TODO: can reduce the range of these loops.
        for (int ist = st + 1; ist < en - 1; ++ist) {
          for (int ien = ist + 1; ien < en; ++ien) {
            UPDATE_CACHE(P, sz, st, energy::TwoLoop(st, en, ist, ien) + arr[ien - ist + 1][ist][P]);
          }
        }
        // Hairpin loops.
        UPDATE_CACHE(P, sz, st, energy::HairpinEnergy(st, en));

        // Multiloops. Look at range [st + 1, en - 1].
        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        auto base_branch_cost = energy::AuGuPenalty(st, en) + multiloop_hack_a + multiloop_hack_b;
        for (int pivsz = HAIRPIN_MIN_SZ + 2; pivsz < sz - 1 - HAIRPIN_MIN_SZ; ++pivsz) {
          UPDATE_CACHE(
              P, sz, st,
              base_branch_cost + arr[pivsz][st + 1][U] +
              arr[sz - pivsz - 2][st + 1 + pivsz][U]);


          // Paired coaxial stacking cases:
          // Left branch ending base, Left branch ending adjacent unpaired base.
          base_t lben = r[st + pivsz - 1], lbenu = r[st + pivsz];

          // st1en1 = leave 1 unpaired at st + 1 and one unpaired at en - 1.
          auto st1en1_cost =
              base_branch_cost + arr[pivsz - 1][st + 2][P] + arr[sz - pivsz - 3][st + 1 + pivsz][U];
          // (.(   )   .) Left left coax - P
          //UPDATE_CACHE(P, sz, st, st1en1_cost + energy::MismatchMediatedCoaxialEnergy());
          // (.(   )3  .) 3' + Left left coax
          // (.(   ).   ) Left right coax
          // (.(   ).  5) 5' + Left right coax
          // (   .(   ).) Left left coax
          // (3  .(   ).) 3' + Left left coax
          // (.   (   ).) Left right coax
          // (.  5(   ).) 5' + Left right coax
          // ((   )   ) Left flush coax
          // ((   )3  ) 3' + Left flush coax
          // ((   )  5) 5' + Left flush coax
          // (   (   )) Right flush coax
          // (  5(   )) 5' + Right flush coax
          // (3  (   )) 3' + Right flush coax
          // (3   ) 3'
          // (   5) 5'
          // (.   .) Terminal mismatch
        }
      }

      // Update unpaired.
      // Choose |st| to be unpaired.
      if (sz)
        UPDATE_CACHE(U, sz, st, arr[sz - 1][st + 1][U]);
      // Pair here.
      auto augu_penalty = energy::AuGuPenalty(st, en);
      UPDATE_CACHE(U, sz, st, arr[sz][st][P] + multiloop_hack_b + augu_penalty);
      // TODO: can reduce range, probably
      for (int pivsz = HAIRPIN_MIN_SZ + 2; pivsz < sz; ++pivsz) {
        augu_penalty = energy::AuGuPenalty(st, st + pivsz - 1);
        // Split.
        UPDATE_CACHE(
            U, sz, st, arr[pivsz][st][P] + multiloop_hack_b +
                       arr[sz - pivsz][st + pivsz][U] + augu_penalty);
        // Choose rest to be empty.
        UPDATE_CACHE(U, sz, st, arr[pivsz][st][P] + multiloop_hack_b + augu_penalty);
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
      energy_t augu_penalty = energy::AuGuPenalty(st, en - 1);
      energy_t value = exterior[st] + arr[en - st][st][P] + augu_penalty;
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
    if (nst == -1) {
      --ext_st;
      continue;
    }
    q.emplace(ext_st - nst, nst, true);
    ext_st = nst;
  }

  while (!q.empty()) {
    int sz, st;
    bool is_paired;
    std::tie(sz, st, is_paired) = q.front();
    int en = st + sz - 1;
    //printf("%d %d %d: %de\n", st, en, int(is_paired), is_paired ? paired[sz][st] : unpaired[sz][st]);
    q.pop();
    if (is_paired) {
      // It's paired, so add it to the folding.
      p[st] = en;
      p[en] = st;

      // Following largely matches the above DP so look up there for comments.
      for (int ist = st + 1; ist < en - 1; ++ist) {
        for (int ien = ist + 1; ien < en; ++ien) {
          int isz = ien - ist + 1;
          auto val = energy::TwoLoop(st, en, ist, ien) + arr[isz][ist][P];
          if (val == arr[sz][st][P]) {
            q.emplace(isz, ist, true);
            goto loopend;
          }
        }
      }

      auto augu_penalty = energy::AuGuPenalty(st, en);
      for (int pivsz = HAIRPIN_MIN_SZ + 2; pivsz < sz - 1; ++pivsz) {
        auto other_pivsz = sz - pivsz - 2;
        auto val = arr[pivsz][st + 1][U] + multiloop_hack_a + multiloop_hack_b +
                   arr[other_pivsz][st + 1 + pivsz][U] + augu_penalty;
        if (val == arr[sz][st][P]) {
          q.emplace(pivsz, st + 1, false);
          q.emplace(other_pivsz, st + 1 + pivsz, false);
          goto loopend;
        }
      }
    } else {
      if (sz && arr[sz - 1][st + 1][U] == arr[sz][st][U]) {
        q.emplace(sz - 1, st + 1, false);
        goto loopend;
      }
      // Pair here.
      auto augu_penalty = energy::AuGuPenalty(st, en);
      if (arr[sz][st][P] + multiloop_hack_b + augu_penalty == arr[sz][st][U]) {
        q.emplace(sz, st, true);
        goto loopend;
      }
      for (int pivsz = HAIRPIN_MIN_SZ + 2; pivsz < sz; ++pivsz) {
        int rem = sz - pivsz;
        augu_penalty = energy::AuGuPenalty(st, st + pivsz - 1);
        if (arr[pivsz][st][P] + multiloop_hack_b +
            arr[rem][st + pivsz][U] + augu_penalty == arr[sz][st][U]) {
          q.emplace(pivsz, st, true);
          q.emplace(rem, st + pivsz, false);
          goto loopend;
        }
        if (arr[pivsz][st][P] + multiloop_hack_b + augu_penalty == arr[sz][st][U]) {
          q.emplace(pivsz, st, true);
          goto loopend;
        }
      }
    }

    loopend:;
  }
  return exterior[N];
}

#undef UPDATE_CACHE

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

energy_t FoldBruteForce(std::unique_ptr<structure::Structure>* s) {
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
  return energy::ComputeEnergy(s);
}

}
}
