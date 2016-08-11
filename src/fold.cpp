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

enum {
  EXT,
  EXT_WC,
  EXT_GU,
  EXT_RCOAX,
  EXT_SIZE
};

#define UPDATE_CACHE(a, value) \
  do { \
    energy_t macro_upd_value_ = (value); \
    if (macro_upd_value_ < CAP_E && macro_upd_value_ < arr[sz][st][a]) { \
      /*printf("Upd %d %d %d %d => %d\n", sz, st, a, arr[sz][st][a], macro_upd_value_);*/ \
      arr[sz][st][a] = macro_upd_value_; \
    } \
  } while (0)

#define UPDATE_EXT(a, na, psz, pst, value) \
  do { \
    energy_t macro_upd_value_ = (value) + exterior[en + 1][na]; \
    if (macro_upd_value_ < CAP_E && macro_upd_value_ < exterior[st][a]) { \
      printf("Ext %d %d %d %d => %d " #value "\n", st, en, a, exterior[st][a], macro_upd_value_); \
      exterior[st][a] = macro_upd_value_; \
      exterior_sts[st][a] = std::make_tuple(psz, pst, en + 1, na); \
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

          // No stacking case. Not canonical TODO switch to U_ST?
          UPDATE_CACHE(P, base_branch_cost + arr[lpivsz][st + 1][U] + arr[rpivsz][st + 1 + lpivsz][U]);
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
      for (int lpivsz = HAIRPIN_MIN_SZ + 2; lpivsz <= sz; ++lpivsz) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        int rpivsz = sz - lpivsz;
        auto pb = r[st + lpivsz - 1], pl1b = r[st + lpivsz - 2];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        auto base00 = arr[lpivsz][st][P] + energy::AuGuPenalty(st, st + lpivsz - 1) + multiloop_hack_b;
        auto base01 = arr[lpivsz - 1][st][P] + energy::AuGuPenalty(st, st + lpivsz - 2) + multiloop_hack_b;
        auto base10 = arr[lpivsz - 1][st + 1][P] + energy::AuGuPenalty(st + 1, st + lpivsz - 1) + multiloop_hack_b;
        auto base11 = arr[lpivsz - 2][st + 1][P] + energy::AuGuPenalty(st + 1, st + lpivsz - 2) + multiloop_hack_b;
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

        // (   ).<(   ). > Right coax forward and backward
        UPDATE_CACHE(U, base01 + arr[rpivsz][st + lpivsz][U_RCOAX]);
        if (st > 0)
          UPDATE_CACHE(U_RCOAX, base01 + energy::MismatchMediatedCoaxialEnergy(
              pl1b, pb, r[st - 1], stb) + right_unpaired);

        // There has to be remaining bases to even have a chance at these cases.
        if (rpivsz > 0) {
          auto pr1b = r[st + lpivsz];
          // (   )<(   ) > Flush coax - U
          UPDATE_CACHE(U, base00 + stacking_e[pb][pr1b][pr1b ^ 2][stb] + arr[rpivsz][st + lpivsz][U_WC]);
          if (pr1b == G || pr1b == U)
            UPDATE_CACHE(U, base00 + stacking_e[pb][pr1b][pr1b ^ 1][stb] + arr[rpivsz][st + lpivsz][U_GU]);
        }
      }

      en++;
    }
  }

  // Exterior loop calculation. There can be no paired base on exterior[en].
  array2d_t<energy_t, EXT_SIZE> exterior(std::size_t(N + 1));
  exterior[N][EXT] = 0;
  // Holds: sz, st, nst, na - nst and na for next index into itself; sz, st for the paired.
  array2d_t<std::tuple<int, int, int, int>, EXT_SIZE> exterior_sts(std::size_t(N), 0xFF);
  for (int st = N - 1; st >= 0; --st) {
    // Case: No pair starting here
    exterior[st][EXT] = exterior[st + 1][EXT];
    for (int sz = HAIRPIN_MIN_SZ + 2; sz < N - st + 1; ++sz) {
      // .   .   .   (   .   .   .   )   <   >
      //           stb  st1b   en1b  enb   rem
      int en = st + sz - 1, rem = N - st - sz;
      auto stb = r[st], st1b = r[st + 1], enb = r[en], en1b = r[en - 1];

      auto base00 = arr[sz][st][P] + energy::AuGuPenalty(st, en);
      auto base01 = arr[sz - 1][st][P] + energy::AuGuPenalty(st, en - 1);
      auto base10 = arr[sz - 1][st + 1][P] + energy::AuGuPenalty(st + 1, en);
      auto base11 = arr[sz - 2][st + 1][P] + energy::AuGuPenalty(st + 1, en - 1);

      // (   )<   >
      UPDATE_EXT(EXT, EXT, sz, st, base00);
      if (IsGu(stb, enb))
        UPDATE_EXT(EXT_GU, EXT, sz, st, base00);
      else
        UPDATE_EXT(EXT_WC, EXT, sz, st, base00);

      // (   )3<   > 3'
      UPDATE_EXT(EXT, EXT, sz - 1, st, base01 + dangle3_e[en1b][enb][stb]);
      // 5(   )<   > 5'
      UPDATE_EXT(EXT, EXT, sz - 1, st + 1, base10 + dangle5_e[enb][stb][st1b]);
      // .(   ).<   > Terminal mismatch
      UPDATE_EXT(EXT, EXT, sz - 2, st + 1, base11 + terminal_e[en1b][enb][stb][st1b]);
      // .(   ).<(   ) > Left coax
      UPDATE_EXT(EXT, EXT, sz - 2, st + 1, base11 + terminal_e[en1b][enb][stb][st1b]);

      // (   ).<(   ). > Right coax forward and backward
      UPDATE_EXT(EXT, EXT_RCOAX, sz - 1, st, base01);
      if (st > 0)
        UPDATE_EXT(EXT_RCOAX, EXT, sz - 1, st, base01 + energy::MismatchMediatedCoaxialEnergy(
            en1b, enb, r[st - 1], stb));

      if (rem > 0) {
        // (   )<(   ) > Flush coax
        auto enrb = r[en + 1];
        UPDATE_EXT(EXT, EXT_WC, sz, st, base00 + stacking_e[enb][enrb][enrb ^ 2][stb]);
        if (enrb == G || enrb == U)
          UPDATE_EXT(EXT, EXT_GU, sz, st, base00 + stacking_e[enb][enrb][enrb ^ 1][stb]);
      }
    }
  }

  // Backtrace.
  p = std::vector<int>(r.size(), -1);
  // sz, st, paired
  std::queue<std::tuple<int, int, int>> q;
  int ext_st = 0, ext_a = EXT;
  while (ext_st < N) {
    int psz, pst, nst, na;
    std::tie(psz, pst, nst, na) = exterior_sts[ext_st][ext_a];
    if (nst == -1) {
      ++ext_st;
      continue;
    }
    //printf("Exterior: %d %d\n", psz, pst);
    q.emplace(psz, pst, P);
    ext_st = nst;
    ext_a = na;
  }

  while (!q.empty()) {
    int sz, st, a;
    std::tie(sz, st, a) = q.front();
    int en = st + sz - 1;
    assert(sz >= HAIRPIN_MIN_SZ + 2);
    auto stb = r[st], st1b = r[st + 1], st2b = r[st + 2], enb = r[en], en1b = r[en - 1], en2b = r[en - 2];
    //printf("%d %d %d: %de\n", st, en, a, arr[sz][st][a]);
    q.pop();
    if (a == P) {
      // It's paired, so add it to the folding.
      p[st] = en;
      p[en] = st;

      // Following largely matches the above DP so look up there for comments.
      for (int ist = st + 1; ist < en - 1; ++ist) {
        for (int ien = ist + 1; ien < en; ++ien) {
          int isz = ien - ist + 1;
          auto val = energy::TwoLoop(st, en, ist, ien) + arr[isz][ist][P];
          if (val == arr[sz][st][P]) {
            q.emplace(isz, ist, P);
            goto loopend;
          }
        }
      }

      auto base_branch_cost = energy::AuGuPenalty(st, en) + multiloop_hack_a + multiloop_hack_b;
      for (int lpivsz = HAIRPIN_MIN_SZ + 2; lpivsz < sz - 3 - HAIRPIN_MIN_SZ; ++lpivsz) {
        auto rpivsz = sz - lpivsz - 2;
        base_t pl1b = r[st + lpivsz - 1], plb = r[st + lpivsz], prb = r[st + lpivsz + 1], pr1b = r[st + lpivsz + 2];

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

        auto val = base_branch_cost + arr[lpivsz][st + 1][U] + arr[rpivsz][st + 1 + lpivsz][U];
        if (val == arr[sz][st][P]) {
          q.emplace(lpivsz, st + 1, U);
          q.emplace(rpivsz, st + 1 + lpivsz, U);
          goto loopend;
        }

        // (.(   )   .) Left left coax - P
        auto outer_coax = energy::MismatchMediatedCoaxialEnergy(stb, st1b, en1b, enb);
        if (l_l10r01_cost + outer_coax == arr[sz][st][P]) {
          q.emplace(lpivsz - 1, st + 2, P);
          q.emplace(rpivsz - 1, st + 1 + lpivsz, U);
          goto loopend;
        }
        // (.(   ).   ) Left right coax
        auto l_inner_coax = energy::MismatchMediatedCoaxialEnergy(pl1b, plb, st1b, st2b);
        if (l_l11r00_cost + l_inner_coax == arr[sz][st][P]) {
          q.emplace(lpivsz - 2, st + 2, P);
          q.emplace(rpivsz, st + 1 + lpivsz, U);
          goto loopend;
        }
        // (   .(   ).) Right left coax
        auto r_inner_coax = energy::MismatchMediatedCoaxialEnergy(en2b, en1b, prb, pr1b);
        if (r_l00r11_cost + r_inner_coax == arr[sz][st][P]) {
          q.emplace(lpivsz, st + 1, U);
          q.emplace(rpivsz - 2, st + 2 + lpivsz, P);
          goto loopend;
        }
        // (.   (   ).) Right right coax
        if (r_l10r01_cost + outer_coax == arr[sz][st][P]) {
          q.emplace(lpivsz - 1, st + 2, U);
          q.emplace(rpivsz - 1, st + 1 + lpivsz, P);
          goto loopend;
        }
        // ((   )   ) Left flush coax
        if (l_l00r00_cost + stacking_e[stb][st1b][plb][enb] == arr[sz][st][P]) {
          q.emplace(lpivsz, st + 1, P);
          q.emplace(rpivsz, st + 1 + lpivsz, U);
          goto loopend;
        }
        // (   (   )) Right flush coax
        if (r_l00r00_cost + stacking_e[stb][prb][en1b][enb] == arr[sz][st][P]) {
          q.emplace(lpivsz, st + 1, U);
          q.emplace(rpivsz, st + 1 + lpivsz, P);
          goto loopend;
        }
        // (3(   )   ) 3'
        if (l_l10r00_cost + dangle3_e[stb][st1b][enb] == arr[sz][st][P]) {
          q.emplace(lpivsz - 1, st + 2, P);
          q.emplace(rpivsz, st + 1 + lpivsz, U);
          goto loopend;
        }
        // ((   )   5) 5'
        if (l_l00r01_cost + dangle5_e[stb][en1b][enb] == arr[sz][st][P]) {
          q.emplace(lpivsz, st + 1, P);
          q.emplace(rpivsz - 1, st + 1 + lpivsz, U);
          goto loopend;
        }
        // (.(   )   .) Terminal mismatch
        if (l_l10r01_cost + terminal_e[stb][st1b][en1b][enb] == arr[sz][st][P]) {
          q.emplace(lpivsz - 1, st + 2, P);
          q.emplace(rpivsz - 1, st + 1 + lpivsz, U);
          goto loopend;
        }
      }
    } else {
      // Unpaired.
      if (sz && arr[sz - 1][st + 1][U] == arr[sz][st][U]) {
        q.emplace(sz - 1, st + 1, U);
        goto loopend;
      }
      // Pair here.
      if (arr[sz][st][P] + multiloop_hack_b + energy::AuGuPenalty(st, en) == arr[sz][st][U]) {
        q.emplace(sz, st, P);
        goto loopend;
      }

      // Pair here.
      for (int lpivsz = HAIRPIN_MIN_SZ + 2; lpivsz <= sz; ++lpivsz) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        int rpivsz = sz - lpivsz;
        auto pb = r[st + lpivsz - 1], pl1b = r[st + lpivsz - 2];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        auto base00 = arr[lpivsz][st][P] + energy::AuGuPenalty(st, st + lpivsz - 1) + multiloop_hack_b;
        auto base01 = arr[lpivsz - 1][st][P] + energy::AuGuPenalty(st, st + lpivsz - 2) + multiloop_hack_b;
        auto base10 = arr[lpivsz - 1][st + 1][P] + energy::AuGuPenalty(st + 1, st + lpivsz - 1) + multiloop_hack_b;
        auto base11 = arr[lpivsz - 2][st + 1][P] + energy::AuGuPenalty(st + 1, st + lpivsz - 2) + multiloop_hack_b;
        // Min is for either placing another unpaired or leaving it as nothing.
        auto right_unpaired = std::min(arr[rpivsz][st + lpivsz][U], 0);

        // (   )<   > - U, U_WC?, U_GU?
        if (base00 + right_unpaired == arr[sz][st][a]) {
          q.emplace(lpivsz, st, P);
          assert((a != U_WC || IsWatsonCrick(stb, pb)) && (a != U_GU || IsGu(stb, pb)));
          if (right_unpaired)
            q.emplace(rpivsz, st + lpivsz, U);
          goto loopend;
        }

        // (   ).<( ** ). > Right coax backward
        if (a == U_RCOAX && st > 0 && base01 + energy::MismatchMediatedCoaxialEnergy(
            pl1b, pb, r[st - 1], stb) + right_unpaired == arr[sz][st][U_RCOAX]) {
          q.emplace(lpivsz - 1, st, P);
          if (right_unpaired)
            q.emplace(rpivsz, st + lpivsz, U);
          goto loopend;
        }

        assert(a == U);

        // (   )3<   > 3' - U
        if (base01 + dangle3_e[pl1b][pb][stb] + right_unpaired == arr[sz][st][U]) {
          q.emplace(lpivsz - 1, st, P);
          if (right_unpaired)
            q.emplace(rpivsz, st + lpivsz, U);
          goto loopend;
        }
        // 5(   )<   > 5' - U
        if (base10 + dangle5_e[pb][stb][st1b] + right_unpaired == arr[sz][st][U]) {
          q.emplace(lpivsz - 1, st + 1, P);
          if (right_unpaired)
            q.emplace(rpivsz, st + lpivsz, U);
          goto loopend;
        }
        // .(   ).<   > Terminal mismatch - U
        if (base11 + terminal_e[pl1b][pb][stb][st1b] + right_unpaired == arr[sz][st][U]) {
          q.emplace(lpivsz - 2, st + 1, P);
          if (right_unpaired)
            q.emplace(rpivsz, st + lpivsz, U);
          goto loopend;
        }
        // .(   ).<(   ) > Left coax - U
        if (base11 + terminal_e[pl1b][pb][stb][st1b] + right_unpaired == arr[sz][st][U]) {
          q.emplace(lpivsz - 2, st + 1, P);
          if (right_unpaired)
            q.emplace(rpivsz, st + lpivsz, U);
          goto loopend;
        }

        // (   ).<(   ). > Right coax forward
        if (base01 + arr[rpivsz][st + lpivsz][U_RCOAX] == arr[sz][st][U]) {
          q.emplace(lpivsz - 1, st, P);
          q.emplace(rpivsz, st + lpivsz, U_RCOAX);
          goto loopend;
        }

        // There has to be remaining bases to even have a chance at these cases.
        if (rpivsz > 0) {
          auto pr1b = r[st + lpivsz];
          // (   )<(   ) > Flush coax - U
          if (base00 + stacking_e[pb][pr1b][pr1b ^ 2][stb] + arr[rpivsz][st + lpivsz][U_WC] == arr[sz][st][U]) {
            q.emplace(lpivsz, st, P);
            q.emplace(rpivsz, st + lpivsz, U_WC);
            goto loopend;
          }
          if ((pr1b == G || pr1b == U) &&
              base00 + stacking_e[pb][pr1b][pr1b ^ 1][stb] + arr[rpivsz][st + lpivsz][U_GU] == arr[sz][st][U]) {
            q.emplace(lpivsz, st, P);
            q.emplace(rpivsz, st + lpivsz, U_GU);
          }
        }
      }
    }

    loopend:;
  }
  return exterior[0][EXT];
}

#undef UPDATE_CACHE
#undef UPDATE_EXT

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
