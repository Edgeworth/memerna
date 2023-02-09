// Copyright 2023 Eliot Courtney.
#include "models/t22/trace/trace.h"

#include <algorithm>
#include <compare>
#include <memory>
#include <stack>
#include <variant>

#include "model/base.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "models/t04/mfe/dp.h"
#include "models/t22/mfe/mfe.h"

namespace mrna::md::t22 {

namespace {

using t04::DP_P;
using t04::DP_U;
using t04::DP_U2;
using t04::DP_U_GU;
using t04::DP_U_RC;
using t04::DP_U_WC;
using t04::EXT;
using t04::EXT_GU;
using t04::EXT_RC;
using t04::EXT_WC;

struct TracebackInternal {
  const Primary& r;
  const Model& em;
  int N;
  const t04::DpArray& dp;
  const t04::ExtArray& ext;
  const Array3D<Energy>& penult;
  TraceResult res;
  std::stack<Index> q;

  TracebackInternal(const Primary& r_, const Model::Ptr& em_, const DpState& state_)
      : r(r_), em(*em_), N(static_cast<int>(r_.size())), dp(state_.t04.dp), ext(state_.t04.ext),
        penult(state_.penult), res((Secondary(N)), Ctds(N)) {}

  bool ComputeTwoStack(int st, int en, int nst, int nen, Energy base, Energy target) {
    // Can either continue by adding a multiloop, adding a internal
    // loop, or ending with a hairpin.
    base += em.penultimate_stack[r[st]][r[nst]][r[nen]][r[en]];

    const auto multiloop = dp[nst + 1][nen - 1][DP_U2] + em.AuGuPenalty(r[nst], r[nen]) + base;
    // TODO(0): set paired values
    if (multiloop == target) {
      q.emplace(t04::Index(nst + 1, nen - 1, DP_U2));
      return true;
    }

    const auto hairpin = em.Hairpin(r, nst, nen) + base;
    if (hairpin == target) return true;

    const int max_inter = std::min(TWOLOOP_MAX_SZ, nen - nst - HAIRPIN_MIN_SZ - 3);
    for (int ist = nst + 1; ist < nst + max_inter + 2; ++ist) {
      for (int ien = nen - max_inter + ist - nst - 2; ien < nen; ++ien) {
        if (dp[ist][ien][DP_P] < CAP_E && ist - nst + nen - ien > 3) {
          const auto twoloop = em.TwoLoop(r, nst, nen, ist, ien) + base + dp[ist][ien][DP_P];
          if (twoloop == target) {
            q.emplace(t04::Index(ist, ien, DP_P));
            return true;
          }
        }
      }
    }
    return false;
  }

  TraceResult Compute() {
    q.emplace(t04::Index(0, -1, EXT));
    while (!q.empty()) {
      auto idx = q.top();
      q.pop();

      if (std::holds_alternative<t04::Index>(idx)) {
        auto t04idx = std::get<t04::Index>(idx);
        int st = t04idx.st;
        int en = t04idx.en;
        int a = t04idx.a;
        // fmt::print("t04idx: st={}, en={}, a={}\n", st, en, a);
        if (en == -1) {
          // fmt::print("  ext: {}\n", ext[st][EXT]);
        } else {
          // fmt::print("  dp: {}\n", dp[st][en][a]);
        }
        if (en == -1) {
          // Case: No pair starting here
          if (a == EXT && st + 1 < N && ext[st + 1][EXT] == ext[st][EXT]) {
            q.emplace(t04::Index(st + 1, -1, EXT));
            goto loopend0;
          }
          for (en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
            // .   .   .   (   .   .   .   )   <   >
            //           stb  st1b   en1b  enb   rem
            const auto stb = r[st];
            const auto st1b = r[st + 1];
            const auto enb = r[en];
            const auto en1b = r[en - 1];
            const auto base00 = dp[st][en][DP_P] + em.AuGuPenalty(stb, enb);
            const auto base01 = dp[st][en - 1][DP_P] + em.AuGuPenalty(stb, en1b);
            const auto base10 = dp[st + 1][en][DP_P] + em.AuGuPenalty(st1b, enb);
            const auto base11 = dp[st + 1][en - 1][DP_P] + em.AuGuPenalty(st1b, en1b);

            // (   )<.( * ). > Right coax backward
            if (a == EXT_RC) {
              // Don't set CTDs here since they will have already been set.
              if (base11 + em.MismatchCoaxial(en1b, enb, stb, st1b) + ext[en + 1][EXT] ==
                  ext[st][EXT_RC]) {
                q.emplace(t04::Index(st + 1, en - 1, DP_P));
                q.emplace(t04::Index(en + 1, -1, EXT));
                goto loopend0;
              }
              continue;
            }

            // (   )<   >
            auto val = base00 + ext[en + 1][EXT];
            if (val == ext[st][a] && (a != EXT_WC || IsWcPair(stb, enb)) &&
                (a != EXT_GU || IsGuPair(stb, enb))) {
              // EXT_WC and EXT_GU will have already had their ctds set.
              if (a == EXT) res.ctd[st] = CTD_UNUSED;
              q.emplace(t04::Index(st, en, DP_P));
              q.emplace(t04::Index(en + 1, -1, EXT));
              goto loopend0;
            }

            // Only look at EXT from here on.
            if (a != EXT) continue;

            // (   )3<   > 3'
            if (base01 + em.dangle3[en1b][enb][stb] + ext[en + 1][EXT] == ext[st][EXT]) {
              res.ctd[st] = CTD_3_DANGLE;
              q.emplace(t04::Index(st, en - 1, DP_P));
              q.emplace(t04::Index(en + 1, -1, EXT));
              goto loopend0;
            }
            // 5(   )<   > 5'
            if (base10 + em.dangle5[enb][stb][st1b] + ext[en + 1][EXT] == ext[st][EXT]) {
              res.ctd[st + 1] = CTD_5_DANGLE;
              q.emplace(t04::Index(st + 1, en, DP_P));
              q.emplace(t04::Index(en + 1, -1, EXT));
              goto loopend0;
            }
            // .(   ).<   > Terminal mismatch
            if (base11 + em.terminal[en1b][enb][stb][st1b] + ext[en + 1][EXT] == ext[st][EXT]) {
              res.ctd[st + 1] = CTD_MISMATCH;
              q.emplace(t04::Index(st + 1, en - 1, DP_P));
              q.emplace(t04::Index(en + 1, -1, EXT));
              goto loopend0;
            }

            if (en < N - 1) {
              // .(   ).<(   ) > Left coax  x
              val = base11 + em.MismatchCoaxial(en1b, enb, stb, st1b);
              if (val + ext[en + 1][EXT_WC] == ext[st][EXT]) {
                res.ctd[st + 1] = CTD_LCOAX_WITH_NEXT;
                res.ctd[en + 1] = CTD_LCOAX_WITH_PREV;
                q.emplace(t04::Index(st + 1, en - 1, DP_P));
                q.emplace(t04::Index(en + 1, -1, EXT_WC));
                goto loopend0;
              }
              if (val + ext[en + 1][EXT_GU] == ext[st][EXT]) {
                res.ctd[st + 1] = CTD_LCOAX_WITH_NEXT;
                res.ctd[en + 1] = CTD_LCOAX_WITH_PREV;
                q.emplace(t04::Index(st + 1, en - 1, DP_P));
                q.emplace(t04::Index(en + 1, -1, EXT_GU));
                goto loopend0;
              }

              // (   )<.(   ). > Right coax forward
              if (base00 + ext[en + 1][EXT_RC] == ext[st][EXT]) {
                res.ctd[st] = CTD_RC_WITH_NEXT;
                res.ctd[en + 2] = CTD_RC_WITH_PREV;
                q.emplace(t04::Index(st, en, DP_P));
                q.emplace(t04::Index(en + 1, -1, EXT_RC));
                goto loopend0;
              }

              // (   )(<   ) > Flush coax
              if (base01 + em.stack[en1b][enb][WcPair(enb)][stb] + ext[en][EXT_WC] ==
                  ext[st][EXT]) {
                res.ctd[st] = CTD_FCOAX_WITH_NEXT;
                res.ctd[en] = CTD_FCOAX_WITH_PREV;
                q.emplace(t04::Index(st, en - 1, DP_P));
                q.emplace(t04::Index(en, -1, EXT_WC));
                goto loopend0;
              }
              if (IsGu(enb) &&
                  base01 + em.stack[en1b][enb][GuPair(enb)][stb] + ext[en][EXT_GU] ==
                      ext[st][EXT]) {
                res.ctd[st] = CTD_FCOAX_WITH_NEXT;
                res.ctd[en] = CTD_FCOAX_WITH_PREV;
                q.emplace(t04::Index(st, en - 1, DP_P));
                q.emplace(t04::Index(en, -1, EXT_GU));
                goto loopend0;
              }
            }
          }
        } else {
          const auto stb = r[st];
          const auto st1b = r[st + 1];
          const auto st2b = r[st + 2];
          const auto enb = r[en];
          const auto en1b = r[en - 1];
          const auto en2b = r[en - 2];
          if (a == DP_P) {
            // It's paired, so add it to the folding.
            res.s[st] = en;
            res.s[en] = st;

            {
              const int max_stack = en - st - HAIRPIN_MIN_SZ + 1;
              const Energy bulge_left = em.Bulge(r, st, en, st + 2, en - 1);
              const Energy bulge_right = em.Bulge(r, st, en, st + 1, en - 2);
              const Energy augu_penalty = em.AuGuPenalty(stb, enb);

              for (int length = 2; 2 * length <= max_stack; ++length) {
                // fmt::print("length = {}, {}\n", length, dp[st][en][DP_P]);
                auto none = penult[st + 1][en - 1][length - 1] +
                    em.stack[r[st]][r[st + 1]][r[en - 1]][r[en]] +
                    em.penultimate_stack[en1b][enb][stb][st1b] + augu_penalty;

                // If we end the stack here, we are done.
                if (length == 2 &&
                    ComputeTwoStack(st, en, st + 1, en - 1, none, dp[st][en][DP_P])) {
                  goto loopend0;
                } else if (none == dp[st][en][DP_P]) {
                  q.emplace(PenultimateIndex(st + 1, en - 1, length - 1));
                  goto loopend0;
                }

                auto left = penult[st + 2][en - 1][length - 1] + bulge_left +
                    em.penultimate_stack[en1b][enb][stb][st2b] + augu_penalty;
                if (length == 2 &&
                    ComputeTwoStack(st, en, st + 2, en - 1, left, dp[st][en][DP_P])) {
                  goto loopend0;
                } else if (left == dp[st][en][DP_P]) {
                  q.emplace(PenultimateIndex(st + 2, en - 1, length - 1));
                  goto loopend0;
                }

                auto right = penult[st + 1][en - 2][length - 1] + bulge_right +
                    em.penultimate_stack[en2b][enb][stb][st1b] + augu_penalty;
                if (length == 2 &&
                    ComputeTwoStack(st, en, st + 1, en - 2, right, dp[st][en][DP_P])) {
                  goto loopend0;
                } else if (right == dp[st][en][DP_P]) {
                  q.emplace(PenultimateIndex(st + 1, en - 2, length - 1));
                  goto loopend0;
                }
              }
            }

            // Following largely matches the above DP so look up there for comments.
            const int max_inter = std::min(TWOLOOP_MAX_SZ, en - st - HAIRPIN_MIN_SZ - 3);
            for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
              for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
                // Try all internal loops. We don't check stacks or 1 nuc bulge loops.
                if (dp[ist][ien][DP_P] < CAP_E && ist - st + en - ien > 3) {
                  const auto val = em.TwoLoop(r, st, en, ist, ien) + dp[ist][ien][DP_P];
                  if (val == dp[st][en][DP_P]) {
                    q.emplace(t04::Index(ist, ien, DP_P));
                    goto loopend0;
                  }
                }
              }
            }

            const auto base_branch_cost =
                em.AuGuPenalty(stb, enb) + em.multiloop_hack_a + em.multiloop_hack_b;
            // (<   ><    >)
            if (base_branch_cost + dp[st + 1][en - 1][DP_U2] == dp[st][en][DP_P]) {
              res.ctd[en] = CTD_UNUSED;
              q.emplace(t04::Index(st + 1, en - 1, DP_U2));
              goto loopend0;
            }
            // (3<   ><   >) 3'
            if (base_branch_cost + dp[st + 2][en - 1][DP_U2] + em.dangle3[stb][st1b][enb] ==
                dp[st][en][DP_P]) {
              res.ctd[en] = CTD_3_DANGLE;
              q.emplace(t04::Index(st + 2, en - 1, DP_U2));
              goto loopend0;
            }
            // (<   ><   >5) 5'
            if (base_branch_cost + dp[st + 1][en - 2][DP_U2] + em.dangle5[stb][en1b][enb] ==
                dp[st][en][DP_P]) {
              res.ctd[en] = CTD_5_DANGLE;
              q.emplace(t04::Index(st + 1, en - 2, DP_U2));
              goto loopend0;
            }
            // (.<   ><   >.) Terminal mismatch
            if (base_branch_cost + dp[st + 2][en - 2][DP_U2] + em.terminal[stb][st1b][en1b][enb] ==
                dp[st][en][DP_P]) {
              res.ctd[en] = CTD_MISMATCH;
              q.emplace(t04::Index(st + 2, en - 2, DP_U2));
              goto loopend0;
            }

            for (int piv = st + HAIRPIN_MIN_SZ + 2; piv < en - HAIRPIN_MIN_SZ - 2; ++piv) {
              const Base pl1b = r[piv - 1];
              const Base plb = r[piv];
              const Base prb = r[piv + 1];
              const Base pr1b = r[piv + 2];

              // (.(   )   .) Left outer coax - P
              const auto outer_coax = em.MismatchCoaxial(stb, st1b, en1b, enb);
              if (base_branch_cost + dp[st + 2][piv][DP_P] + em.multiloop_hack_b +
                      em.AuGuPenalty(st2b, plb) + dp[piv + 1][en - 2][DP_U] + outer_coax ==
                  dp[st][en][DP_P]) {
                res.ctd[en] = CTD_LCOAX_WITH_NEXT;
                res.ctd[st + 2] = CTD_LCOAX_WITH_PREV;
                q.emplace(t04::Index(st + 2, piv, DP_P));
                q.emplace(t04::Index(piv + 1, en - 2, DP_U));
                goto loopend0;
              }
              // (.   (   ).) Right outer coax
              if (base_branch_cost + dp[st + 2][piv][DP_U] + em.multiloop_hack_b +
                      em.AuGuPenalty(prb, en2b) + dp[piv + 1][en - 2][DP_P] + outer_coax ==
                  dp[st][en][DP_P]) {
                res.ctd[en] = CTD_RC_WITH_PREV;
                res.ctd[piv + 1] = CTD_RC_WITH_NEXT;
                q.emplace(t04::Index(st + 2, piv, DP_U));
                q.emplace(t04::Index(piv + 1, en - 2, DP_P));
                goto loopend0;
              }

              // (.(   ).   ) Left inner coax
              if (base_branch_cost + dp[st + 2][piv - 1][DP_P] + em.multiloop_hack_b +
                      em.AuGuPenalty(st2b, pl1b) + dp[piv + 1][en - 1][DP_U] +
                      em.MismatchCoaxial(pl1b, plb, st1b, st2b) ==
                  dp[st][en][DP_P]) {
                res.ctd[en] = CTD_RC_WITH_NEXT;
                res.ctd[st + 2] = CTD_RC_WITH_PREV;
                q.emplace(t04::Index(st + 2, piv - 1, DP_P));
                q.emplace(t04::Index(piv + 1, en - 1, DP_U));
                goto loopend0;
              }
              // (   .(   ).) Right inner coax
              if (base_branch_cost + dp[st + 1][piv][DP_U] + em.multiloop_hack_b +
                      em.AuGuPenalty(pr1b, en2b) + dp[piv + 2][en - 2][DP_P] +
                      em.MismatchCoaxial(en2b, en1b, prb, pr1b) ==
                  dp[st][en][DP_P]) {
                res.ctd[en] = CTD_LCOAX_WITH_PREV;
                res.ctd[piv + 2] = CTD_LCOAX_WITH_NEXT;
                q.emplace(t04::Index(st + 1, piv, DP_U));
                q.emplace(t04::Index(piv + 2, en - 2, DP_P));
                goto loopend0;
              }

              // ((   )   ) Left flush coax
              if (base_branch_cost + dp[st + 1][piv][DP_P] + em.multiloop_hack_b +
                      em.AuGuPenalty(st1b, plb) + dp[piv + 1][en - 1][DP_U] +
                      em.stack[stb][st1b][plb][enb] ==
                  dp[st][en][DP_P]) {
                res.ctd[en] = CTD_FCOAX_WITH_NEXT;
                res.ctd[st + 1] = CTD_FCOAX_WITH_PREV;
                q.emplace(t04::Index(st + 1, piv, DP_P));
                q.emplace(t04::Index(piv + 1, en - 1, DP_U));
                goto loopend0;
              }
              // (   (   )) Right flush coax
              if (base_branch_cost + dp[st + 1][piv][DP_U] + em.multiloop_hack_b +
                      em.AuGuPenalty(prb, en1b) + dp[piv + 1][en - 1][DP_P] +
                      em.stack[stb][prb][en1b][enb] ==
                  dp[st][en][DP_P]) {
                res.ctd[en] = CTD_FCOAX_WITH_PREV;
                res.ctd[piv + 1] = CTD_FCOAX_WITH_NEXT;
                q.emplace(t04::Index(st + 1, piv, DP_U));
                q.emplace(t04::Index(piv + 1, en - 1, DP_P));
                goto loopend0;
              }
            }

            // Done with paired. We might not have jumped to loopend if this was a hairpin.
            continue;
          }

          // Deal with the rest of the cases:
          // Left unpaired. Either DP_U or DP_U2.
          if (st + 1 < en && (a == DP_U || a == DP_U2) && dp[st + 1][en][a] == dp[st][en][a]) {
            q.emplace(t04::Index(st + 1, en, a));
            goto loopend0;
          }

          // Pair here.
          for (int piv = st + HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
            //   (   .   )<   (
            // stb pl1b pb   pr1b
            const auto pb = r[piv];
            const auto pl1b = r[piv - 1];
            // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the
            // right.
            const auto base00 = dp[st][piv][DP_P] + em.AuGuPenalty(stb, pb) + em.multiloop_hack_b;
            const auto base01 =
                dp[st][piv - 1][DP_P] + em.AuGuPenalty(stb, pl1b) + em.multiloop_hack_b;
            const auto base10 =
                dp[st + 1][piv][DP_P] + em.AuGuPenalty(st1b, pb) + em.multiloop_hack_b;
            const auto base11 =
                dp[st + 1][piv - 1][DP_P] + em.AuGuPenalty(st1b, pl1b) + em.multiloop_hack_b;

            // Min is for either placing another unpaired or leaving it as nothing.
            // If we're at U2, don't allow leaving as nothing.
            auto right_unpaired = dp[piv + 1][en][DP_U];
            if (a != DP_U2) right_unpaired = std::min(right_unpaired, ZERO_E);

            // Check a == U_RC:
            // (   )<.( ** ). > Right coax backward
            if (a == DP_U_RC) {
              if (base11 + em.MismatchCoaxial(pl1b, pb, stb, st1b) + right_unpaired ==
                  dp[st][en][DP_U_RC]) {
                // Ctds were already set from the recurrence that called this.
                q.emplace(t04::Index(st + 1, piv - 1, DP_P));
                if (right_unpaired != ZERO_E) q.emplace(t04::Index(piv + 1, en, DP_U));
                goto loopend0;
              }
              continue;
            }

            // (   )<   > - U, U2, U_WC?, U_GU?
            if (base00 + right_unpaired == dp[st][en][a] && (a != DP_U_WC || IsWcPair(stb, pb)) &&
                (a != DP_U_GU || IsGuPair(stb, pb))) {
              // If U_WC, or U_GU, we were involved in some sort of coaxial stack previously, and
              // were already set.
              if (a != DP_U_WC && a != DP_U_GU) res.ctd[st] = CTD_UNUSED;
              q.emplace(t04::Index(st, piv, DP_P));
              if (a == DP_U2 || right_unpaired != ZERO_E) q.emplace(t04::Index(piv + 1, en, DP_U));
              goto loopend0;
            }

            // The rest of the cases are for U and U2.
            if (a != DP_U && a != DP_U2) continue;

            // (   )3<   > 3' - U, U2
            if (base01 + em.dangle3[pl1b][pb][stb] + right_unpaired == dp[st][en][a]) {
              res.ctd[st] = CTD_3_DANGLE;
              q.emplace(t04::Index(st, piv - 1, DP_P));
              if (a == DP_U2 || right_unpaired != ZERO_E) q.emplace(t04::Index(piv + 1, en, DP_U));
              goto loopend0;
            }
            // 5(   )<   > 5' - U, U2
            if (base10 + em.dangle5[pb][stb][st1b] + right_unpaired == dp[st][en][a]) {
              res.ctd[st + 1] = CTD_5_DANGLE;
              q.emplace(t04::Index(st + 1, piv, DP_P));
              if (a == DP_U2 || right_unpaired != ZERO_E) q.emplace(t04::Index(piv + 1, en, DP_U));
              goto loopend0;
            }
            // .(   ).<   > Terminal mismatch - U, U2
            if (base11 + em.terminal[pl1b][pb][stb][st1b] + right_unpaired == dp[st][en][a]) {
              res.ctd[st + 1] = CTD_MISMATCH;
              q.emplace(t04::Index(st + 1, piv - 1, DP_P));
              if (a == DP_U2 || right_unpaired != ZERO_E) q.emplace(t04::Index(piv + 1, en, DP_U));
              goto loopend0;
            }
            // .(   ).<(   ) > Left coax - U, U2
            auto val = base11 + em.MismatchCoaxial(pl1b, pb, stb, st1b);
            if (val + dp[piv + 1][en][DP_U_WC] == dp[st][en][a]) {
              res.ctd[st + 1] = CTD_LCOAX_WITH_NEXT;
              res.ctd[piv + 1] = CTD_LCOAX_WITH_PREV;
              q.emplace(t04::Index(st + 1, piv - 1, DP_P));
              q.emplace(t04::Index(piv + 1, en, DP_U_WC));
              goto loopend0;
            }
            if (val + dp[piv + 1][en][DP_U_GU] == dp[st][en][a]) {
              res.ctd[st + 1] = CTD_LCOAX_WITH_NEXT;
              res.ctd[piv + 1] = CTD_LCOAX_WITH_PREV;
              q.emplace(t04::Index(st + 1, piv - 1, DP_P));
              q.emplace(t04::Index(piv + 1, en, DP_U_GU));
              goto loopend0;
            }

            // (   )<.(   ). > Right coax forward - U, U2
            if (base00 + dp[piv + 1][en][DP_U_RC] == dp[st][en][a]) {
              res.ctd[st] = CTD_RC_WITH_NEXT;
              res.ctd[piv + 2] = CTD_RC_WITH_PREV;
              q.emplace(t04::Index(st, piv, DP_P));
              q.emplace(t04::Index(piv + 1, en, DP_U_RC));
              goto loopend0;
            }

            // (   )(<   ) > Flush coax - U, U2
            if (base01 + em.stack[pl1b][pb][WcPair(pb)][stb] + dp[piv][en][DP_U_WC] ==
                dp[st][en][a]) {
              res.ctd[st] = CTD_FCOAX_WITH_NEXT;
              res.ctd[piv] = CTD_FCOAX_WITH_PREV;
              q.emplace(t04::Index(st, piv - 1, DP_P));
              q.emplace(t04::Index(piv, en, DP_U_WC));
              goto loopend0;
            }
            if ((IsGu(pb)) &&
                base01 + em.stack[pl1b][pb][GuPair(pb)][stb] + dp[piv][en][DP_U_GU] ==
                    dp[st][en][a]) {
              res.ctd[st] = CTD_FCOAX_WITH_NEXT;
              res.ctd[piv] = CTD_FCOAX_WITH_PREV;
              q.emplace(t04::Index(st, piv - 1, DP_P));
              q.emplace(t04::Index(piv, en, DP_U_GU));
              goto loopend0;
            }
          }
        }
      loopend0 : {};
      } else {
        auto penult_idx = std::get<PenultimateIndex>(idx);
        const int st = penult_idx.st;
        const int en = penult_idx.en;
        const int length = penult_idx.len;
        Energy none = ZERO_E;
        Energy left = ZERO_E;
        Energy right = ZERO_E;

        // It's paired, so add it to the folding.
        res.s[st] = en;
        res.s[en] = st;
        // fmt::print("penult: {} {} {} {}\n", st, en, length, penult[st][en][length]);

        const auto bulge_left = em.Bulge(r, st, en, st + 2, en - 1);
        const auto bulge_right = em.Bulge(r, st, en, st + 1, en - 2);

        none = penult[st + 1][en - 1][length - 1] + em.stack[r[st]][r[st + 1]][r[en - 1]][r[en]];
        if (length == 2 && ComputeTwoStack(st, en, st + 1, en - 1, none, penult[st][en][length])) {
          goto loopend1;
        } else if (none == penult[st][en][length]) {
          q.emplace(PenultimateIndex(st + 1, en - 1, length - 1));
          goto loopend1;
        }

        left = penult[st + 2][en - 1][length - 1] + bulge_left;
        if (length == 2 && ComputeTwoStack(st, en, st + 2, en - 1, left, penult[st][en][length])) {
          goto loopend1;
        } else if (left == penult[st][en][length]) {
          q.emplace(PenultimateIndex(st + 2, en - 1, length - 1));
          goto loopend1;
        }
        right = penult[st + 1][en - 2][length - 1] + bulge_right;
        if (length == 2 && ComputeTwoStack(st, en, st + 1, en - 2, right, penult[st][en][length])) {
          goto loopend1;
        } else if (right == penult[st][en][length]) {
          q.emplace(PenultimateIndex(st + 1, en - 2, length - 1));
          goto loopend1;
        }

        // Should not reach here.
        bug();

      loopend1 : {};
      }
    }
    return std::move(res);
  }
};

}  // namespace

TraceResult Traceback(const Primary& r, const Model::Ptr& em, const DpState& state) {
  return TracebackInternal(r, em, state).Compute();
}

}  // namespace mrna::md::t22
