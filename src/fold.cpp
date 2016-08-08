#include <queue>
#include "fold.h"

namespace memerna {
namespace fold {

enum {
  P,  // For the paired array.
  U,  // For the unpaired array.
  U_WC,  // Unpaired but must start with a branch not involved in a CTD interaction that is not GU.
  U_GU,  // Unpaired but must start with a branch not involved in a CTD interaction that is GU.
  U_RCOAX,  // Unpaired but must start with a branch involved in a right coaxial stack - includes energy for it.
  ARR_SIZE
};

#define UPDATE_CACHE(a, value) \
  do { \
    energy_t macro_upd_value_ = (value); \
    if (macro_upd_value_ < CAP_E && macro_upd_value_ < arr[sz][st][a]) { \
      arr[sz][st][a] = macro_upd_value_; \
    } \
  } while (0)


energy_t Fold(std::unique_ptr<structure::Structure>* s) {
  int N = int(r.size());
  assert(N > 0);
  // Automatically initialised to MAX_E.
  array3d_t<energy_t, ARR_SIZE> arr(r.size() + 1);
  assert(arr[N - 1][0][U_WC] == MAX_E);

  // TODO: check bounds
  // TODO: au/gu penalty?
  // sz includes st and en.
  static_assert(HAIRPIN_MIN_SZ >= 2, "Minimum hairpin size >= 2 is relied upon in some expressions.");
  for (int sz = HAIRPIN_MIN_SZ + 2; sz <= N; ++sz) {
    int en = sz - 1;
    for (int st = 0; st < N - sz + 1; ++st) {
      base_t stb = r[st], st1b = r[st + 1], st2b = r[st + 2], enb = r[en], en1b = r[en - 1], en2b = r[en - 2];

      // Update paired - only if can actually pair.
      if (CanPair(r[st], r[en])) {
        // Internal loops, bulge loops, and stacking. TODO: Lyngso's
        // TODO: can reduce the range of these loops.
        for (int ist = st + 1; ist < en - 1; ++ist) {
          for (int ien = ist + 1; ien < en; ++ien) {
            UPDATE_CACHE(P, energy::TwoLoop(st, en, ist, ien) + arr[ien - ist + 1][ist][P]);
          }
        }
        // Hairpin loops.
        UPDATE_CACHE(P, energy::HairpinEnergy(st, en));

        // Multiloops. Look at range [st + 1, en - 1].
        // Cost for initiation + one branch. Include AU/GU penalty for ending multiloop helix.
        auto base_branch_cost = energy::AuGuPenalty(st, en) + multiloop_hack_a + multiloop_hack_b;
        for (int lpivsz = HAIRPIN_MIN_SZ + 2; lpivsz < sz - 3 - HAIRPIN_MIN_SZ; ++lpivsz) {
          int rpivsz = sz - lpivsz - 2;
          // No stacking case. Not canonical TODO switch to U_ST?
          UPDATE_CACHE(P, base_branch_cost + arr[lpivsz][st + 1][U] + arr[rpivsz][st + 1 + lpivsz][U]);

          // Paired coaxial stacking cases:
          base_t pl1b = r[st + lpivsz - 1], plb = r[st + lpivsz], prb = r[st + lpivsz + 1], pr1b = r[st + lpivsz + 2];

          //   (   .   (   .   .   .   )   .   |   .   (   .   .   .   )   .   )
          // stb st1b st2b          pl1b  plb     prb  pr1b         en2b en1b enb

          // l_lABrCD =
          //   l => place a paired on the left
          //   skip A bases on left of left branch, B bases on right of left branch
          //   skip C bases on left of right branch, D bases on right of right branch.
          auto l_l00r00_cost =
              base_branch_cost + arr[lpivsz][st + 1][P] + arr[rpivsz][st + 1 + lpivsz][U];
          auto l_l00r01_cost =
              base_branch_cost + arr[lpivsz][st + 1][P] + arr[rpivsz - 1][st + 1 + lpivsz][U];
          auto l_l10r00_cost =
              base_branch_cost + arr[lpivsz - 1][st + 2][P] + arr[rpivsz][st + 1 + lpivsz][U];
          auto l_l10r01_cost =
              base_branch_cost + arr[lpivsz - 1][st + 2][P] + arr[rpivsz - 1][st + 1 + lpivsz][U];
          auto l_l11r00_cost =
              base_branch_cost + arr[lpivsz - 2][st + 2][P] + arr[rpivsz][st + 1 + lpivsz][U];

          auto r_l00r00_cost =
              base_branch_cost + arr[lpivsz][st + 1][U] + arr[rpivsz][st + 1 + lpivsz][P];
          auto r_l10r01_cost =
              base_branch_cost + arr[lpivsz - 1][st + 2][U] + arr[rpivsz - 1][st + 1 + lpivsz][P];
          auto r_l00r11_cost =
              base_branch_cost + arr[lpivsz][st + 1][U] + arr[rpivsz - 2][st + 2 + lpivsz][P];

          // (.(   )   .) Left left coax - P
          auto outer_coax = energy::MismatchMediatedCoaxialEnergy(stb, st1b, en1b, enb);
          UPDATE_CACHE(P, l_l10r01_cost + outer_coax);
          // (.(   ).   ) Left right coax
          auto l_inner_coax = energy::MismatchMediatedCoaxialEnergy(pl1b, plb, st1b, st2b);
          UPDATE_CACHE(P, l_l11r00_cost + l_inner_coax);
          // (   .(   ).) Right left coax
          auto r_inner_coax = energy::MismatchMediatedCoaxialEnergy(en2b, en1b, prb, pr1b);
          UPDATE_CACHE(P, r_l00r11_cost + r_inner_coax);
          // (.   (   ).) Right right coax
          UPDATE_CACHE(P, r_l10r01_cost + outer_coax);
          // ((   )   ) Left flush coax
          UPDATE_CACHE(P, l_l00r00_cost + stacking_e[stb][st1b][plb][enb]);
          // (   (   )) Right flush coax
          UPDATE_CACHE(P, r_l00r00_cost + stacking_e[stb][prb][en1b][enb]);
          // (3(   )   ) 3'
          UPDATE_CACHE(P, l_l10r00_cost + dangle3_e[stb][st1b][enb]);
          // ((   )   5) 5'
          UPDATE_CACHE(P, l_l00r01_cost + dangle5_e[stb][en1b][enb]);
          // (.(   )   .) Terminal mismatch
          UPDATE_CACHE(P, l_l10r01_cost + terminal_e[stb][st1b][en1b][enb]);
        }
      }

      // Update unpaired.
      // Choose |st| to be unpaired.
      if (sz)
        UPDATE_CACHE(U, arr[sz - 1][st + 1][U]);
      // Pair here.
      UPDATE_CACHE(U, arr[sz][st][P] + multiloop_hack_b + energy::AuGuPenalty(st, en));
      for (int lpivsz = HAIRPIN_MIN_SZ + 2; lpivsz < sz; ++lpivsz) {
        int rpivsz = sz - lpivsz;

        //   (   .   )<   (
        // stb pl1b pb   pr1b

        auto pb = r[st + lpivsz - 1], pl1b = r[st + lpivsz - 2], pr1b = r[st + lpivsz];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        auto base00 = arr[lpivsz][st][P] + energy::AuGuPenalty(st, st + lpivsz - 1) + multiloop_hack_b;
        auto base01 = arr[lpivsz - 1][st][P] + energy::AuGuPenalty(st, st + lpivsz - 2) + multiloop_hack_b;
        auto base10 = arr[lpivsz - 1][st + 1][P] + energy::AuGuPenalty(st + 1, st + lpivsz - 2) + multiloop_hack_b;
        auto base11 = arr[lpivsz - 2][st + 1][P] + energy::AuGuPenalty(st + 1, st + lpivsz - 3) + multiloop_hack_b;
        // Min is for either placing another unpaired or leaving it as nothing.
        auto right_unpaired = std::min(arr[rpivsz][st + lpivsz][U], 0);

        // (   )<   > - U, U_WC?, U_GU?
        auto val = base00 + right_unpaired;
        UPDATE_CACHE(U, val);
        if (IsGu(stb, pb))
          UPDATE_CACHE(U_GU, val);
        else
          UPDATE_CACHE(U_WC, val);

        // (   )3<   > 3' - U
        UPDATE_CACHE(U, base01 + dangle3_e[pl1b][pb][stb] + right_unpaired);

        // 5(   )<   > 5' - U
        UPDATE_CACHE(U, base10 + dangle5_e[pb][stb][st1b] + right_unpaired);

        // .(   ).<   > Terminal mismatch - U
        UPDATE_CACHE(U, base11 + terminal_e[pl1b][pb][stb][st1b] + right_unpaired);

        // .(   ).<(   ) > Left coax - U
        UPDATE_CACHE(U, base11 + terminal_e[pl1b][pb][stb][st1b] + right_unpaired);

        // (   )<(   ) > Flush coax - U
        UPDATE_CACHE(U, base00 + stacking_e[pb][pr1b][pr1b ^ 2][stb] + arr[rpivsz][st + lpivsz][U_WC]);
        if (pr1b == G || pr1b == U)
          UPDATE_CACHE(U, base00 + stacking_e[pb][pr1b][pr1b ^ 1][stb] + arr[rpivsz][st + lpivsz][U_GU]);

        // (   ).<(   ). > Right coax forward and backward
        UPDATE_CACHE(U, base01 + arr[rpivsz][st + lpivsz][U_RCOAX]);
        if (st > 0)
          UPDATE_CACHE(U_RCOAX, base01 + energy::MismatchMediatedCoaxialEnergy(
              pl1b, pb, r[st - 1], stb) + right_unpaired);
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
      energy_t penalty = energy::AuGuPenalty(st, en - 1);
      energy_t value = exterior[st] + arr[en - st][st][P] + penalty;
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

      auto penalty = energy::AuGuPenalty(st, en);
      for (int pivsz = HAIRPIN_MIN_SZ + 2; pivsz < sz - 1; ++pivsz) {
        auto other_pivsz = sz - pivsz - 2;
        auto val = arr[pivsz][st + 1][U] + multiloop_hack_a + multiloop_hack_b +
                   arr[other_pivsz][st + 1 + pivsz][U] + penalty;
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
      auto penalty = energy::AuGuPenalty(st, en);
      if (arr[sz][st][P] + multiloop_hack_b + penalty == arr[sz][st][U]) {
        q.emplace(sz, st, true);
        goto loopend;
      }
      for (int pivsz = HAIRPIN_MIN_SZ + 2; pivsz < sz; ++pivsz) {
        int rem = sz - pivsz;
        penalty = energy::AuGuPenalty(st, st + pivsz - 1);
        if (arr[pivsz][st][P] + multiloop_hack_b +
            arr[rem][st + pivsz][U] + penalty == arr[sz][st][U]) {
          q.emplace(pivsz, st, true);
          q.emplace(rem, st + pivsz, false);
          goto loopend;
        }
        if (arr[pivsz][st][P] + multiloop_hack_b + penalty == arr[sz][st][U]) {
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
