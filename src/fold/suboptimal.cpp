#include <set>
#include "suboptimal.h"
#include "parsing.h"

namespace memerna {
namespace fold {

struct suboptimal_node_t {
  // TODO Since nodes form a tree, can save memory on these two vectors.
  // TODO Will also save time in the comparison.
  // State should be fully defined by |not_yet_expanded| and |history|, which denote
  // what it has done so far, and what it can do from now.
  std::vector<std::tuple<int, int, int>> not_yet_expanded;
  std::vector<std::tuple<int, int, int>> history;
  energy_t energy;  // Stores the minimum energy this state could have.

  bool operator<(const suboptimal_node_t& o) const {
    if (energy != o.energy) return energy < o.energy;
    if (not_yet_expanded != o.not_yet_expanded) return not_yet_expanded < o.not_yet_expanded;
    if (history != o.history) return history < o.history;
    assert(false);  // Should never happen.
    return false;
  }
};

#define PRUNE_INSERT(prune_, node_) \
  do { \
    if (node_.energy <= max_energy) { \
      if (int(prune_.size()) >= max_structures && (--prune_.end())->energy > node_.energy) \
        prune_.erase(--prune_.end()); \
      if (int(prune_.size()) < max_structures) \
        prune_.insert(node_); \
    } \
  } while (0)

// Macro scumming.
#define UNPACK(...) __VA_ARGS__
// Creates and inserts a new node with energy |energy| that doesn't
// need to expand any more ranges than it currently has.
#define EXPAND0(energy_) \
  do { \
    newnode.energy = (energy_); \
    PRUNE_INSERT(q, newnode); \
  } while (0)
// Creates and inserts a new node with energy |energy| that needs to expand the given ranges.
#define EXPAND1(energy_, nye_) \
  do { \
    newnode.not_yet_expanded.emplace_back(UNPACK nye_); \
    newnode.energy = (energy_); \
    PRUNE_INSERT(q, newnode); \
    newnode.not_yet_expanded.pop_back(); \
  } while (0)
// Creates and inserts a new node with energy |energy| that needs to expand the two given ranges.
#define EXPAND2(energy_, nye0_, nye1_) \
  do { \
    newnode.not_yet_expanded.emplace_back(UNPACK nye0_); \
    newnode.not_yet_expanded.emplace_back(UNPACK nye1_); \
    newnode.energy = (energy_); \
    PRUNE_INSERT(q, newnode); \
    newnode.not_yet_expanded.pop_back(); \
    newnode.not_yet_expanded.pop_back(); \
  } while (0)

std::vector<computed_t> SuboptimalTraceback0(
    energy_t max_energy, int max_structures,
    const array3d_t<energy_t, DP_SIZE>& arr,
    const array2d_t<energy_t, EXT_SIZE>& exterior) {
  // Basic idea of suboptimal traceback is look at all possible choices from a state, and expand just one of them.
  // Fully expanding one of them means there will be no duplicates in the tree.
  // Cull the ones not inside the window or when we have more than |max_structures|.
  // We don't have to check for expanding impossible states indirectly, since they will have MAX_E, be
  // above max_energy, and be instantly culled (callers use CAP_E for no energy limit).
  verify_expr(max_structures > 0, "must request at least one structure");
  int N = int(r.size());
  std::set<suboptimal_node_t> finished;
  std::set<suboptimal_node_t> q;

  suboptimal_node_t newnode = {{std::make_tuple(0, -1, EXT)}, {}, exterior[0][EXT]};
  q.insert(newnode);
  while (!q.empty()) {
//    printf("q size: %zu\n", q.size());
//    for (const auto& node : q) {
//      printf("  %de %zu: ", node.energy, node.not_yet_expanded.size());
//      for (const auto& nye : node.not_yet_expanded)
//        printf("(%d, %d, %d), ", std::get<0>(nye), std::get<1>(nye), std::get<2>(nye));
//      printf("\n");
//    }
    auto node = std::move(*q.begin());
    q.erase(q.begin());
    // Finished state.
    if (node.not_yet_expanded.empty()) {
      PRUNE_INSERT(finished, node);
      continue;
    }

    // If we found a non-finished node, but |finished| is full, and the worst in |finished| is as good as
    // our current node (which is the best in |q|), then we can exit.
    if (int(finished.size()) >= max_structures && (--finished.end())->energy <= node.energy)
      break;

    const auto& to_expand = node.not_yet_expanded.back();
    node.not_yet_expanded.pop_back();
    node.history.push_back(to_expand);  // Add to history.
    int st, en, a;
    std::tie(st, en, a) = to_expand;

    // Initialise, since we only make small modifications to it, using the EXPAND macros.
    newnode = node;
    // Temporary variable to hold energy calculations.
    energy_t energy = 0;

    // Exterior loop
    if (en == -1) {
      // We try replacing what we do at (st, a) with a bunch of different cases, so we use this energy as a base.
      energy_t base_energy = node.energy - exterior[st][a];
      if (a == EXT) {
        // Base case: do nothing.
        if (st == N)
          EXPAND0(base_energy);
        else
          // Case: No pair starting here (for EXT only)
          EXPAND1(base_energy + exterior[st + 1][EXT], (st + 1, -1, EXT));
      }
      for (en = st + constants::HAIRPIN_MIN_SZ + 1; en < N; ++en) {
        // .   .   .   (   .   .   .   )   <   >
        //           stb  st1b   en1b  enb   rem
        auto stb = r[st], st1b = r[st + 1], enb = r[en], en1b = r[en - 1];

        auto base00 = arr[st][en][DP_P] + energy::AuGuPenalty(stb, enb);
        auto base01 = arr[st][en - 1][DP_P] + energy::AuGuPenalty(stb, en1b);
        auto base10 = arr[st + 1][en][DP_P] + energy::AuGuPenalty(st1b, enb);
        auto base11 = arr[st + 1][en - 1][DP_P] + energy::AuGuPenalty(st1b, en1b);

        newnode = node;

        // (   ).<( * ). > Right coax backward
        if (st > 0 && a == EXT_RCOAX) {
          energy = base_energy + base01 + energy::MismatchCoaxial(
              en1b, enb, r[st - 1], stb) + exterior[en + 1][EXT];
          EXPAND2(energy, (en + 1, -1, EXT), (st, en - 1, DP_P));
        }

        // Cases for EXT, EXT_WC, EXT_GU.
        if (a == EXT || (a == EXT_WC && IsWatsonCrick(stb, enb)) || (a == EXT_GU && IsGu(stb, enb))) {
          // (   )<   >
          energy = base_energy + base00 + exterior[en + 1][EXT];
          EXPAND2(energy, (en + 1, -1, EXT), (st, en, DP_P));
          // (   )3<   > 3'
          energy = base_energy + base01 + g_dangle3[en1b][enb][stb] + exterior[en + 1][EXT];
          EXPAND2(energy, (en + 1, -1, EXT), (st, en - 1, DP_P));
        }

        // Everything after this is only for EXT.
        if (a != EXT) continue;

        // 5(   )<   > 5'
        energy = base_energy + base10 + g_dangle5[enb][stb][st1b] + exterior[en + 1][EXT];
        EXPAND2(energy, (en + 1, -1, EXT), (st + 1, en, DP_P));

        // .(   ).<   > Terminal mismatch
        energy = base_energy + base11 + g_terminal[en1b][enb][stb][st1b] + exterior[en + 1][EXT];
        EXPAND2(energy, (en + 1, -1, EXT), (st + 1, en - 1, DP_P));

        // .(   ).<(   ) > Left coax
        energy = base_energy + base11 + energy::MismatchCoaxial(en1b, enb, stb, st1b);
        EXPAND2(energy + exterior[en + 1][EXT_GU], (en + 1, -1, EXT_GU), (st + 1, en - 1, DP_P));
        EXPAND2(energy + exterior[en + 1][EXT_WC], (en + 1, -1, EXT_WC), (st + 1, en - 1, DP_P));

        // (   ).<(   ). > Right coax forward
        energy = base_energy + base01 + exterior[en + 1][EXT_RCOAX];
        EXPAND2(energy, (en + 1, -1, EXT_RCOAX), (st, en - 1, DP_P));

        if (en < N - 1) {
          // (   )<(   ) > Flush coax
          auto enrb = r[en + 1];
          energy = base_energy + base00 + g_stack[enb][enrb][enrb ^ 3][stb] + exterior[en + 1][EXT_WC];
          EXPAND2(energy, (en + 1, -1, EXT_WC), (st, en, DP_P));

          if (enrb == G || enrb == U) {
            energy = base_energy + base00 + g_stack[enb][enrb][enrb ^ 1][stb] + exterior[en + 1][EXT_GU];
            EXPAND2(energy, (en + 1, -1, EXT_GU), (st, en, DP_P));
          }
        }
      }
      // Finished exterior loop, don't do anymore.
      continue;
    }

    // Subtract the minimum energy of the contribution at this node.
    energy_t base_energy = node.energy - arr[st][en][a];
    // Declare the usual base aliases.
    auto stb = r[st], st1b = r[st + 1], st2b = r[st + 2], enb = r[en], en1b = r[en - 1], en2b = r[en - 2];

    // Normal stuff
    if (a == DP_P) {
      // Two loops.
      int max_inter = std::min(constants::TWOLOOP_MAX_SZ, en - st - constants::HAIRPIN_MIN_SZ - 3);
      for (int ist = st + 1; ist < st + max_inter + 2; ++ist) {
        for (int ien = en - max_inter + ist - st - 2; ien < en; ++ien) {
          energy = base_energy + energy::TwoLoop(st, en, ist, ien) + arr[ist][ien][DP_P];
          EXPAND1(energy, (ist, ien, DP_P));
        }
      }

      // Hairpin loop
      energy = base_energy + energy::Hairpin(st, en);
      EXPAND0(energy);

      auto base_and_branch = base_energy + energy::AuGuPenalty(stb, enb) + g_multiloop_hack_a + g_multiloop_hack_b;
      // (<   ><    >)
      energy = base_and_branch + arr[st + 1][en - 1][DP_U2];
      EXPAND1(energy, (st + 1, en - 1, DP_U2));
      // (3<   ><   >) 3'
      energy = base_and_branch + arr[st + 2][en - 1][DP_U2] + g_dangle3[stb][st1b][enb];
      EXPAND1(energy, (st + 2, en - 1, DP_U2));
      // (<   ><   >5) 5'
      energy = base_and_branch + arr[st + 1][en - 2][DP_U2] + g_dangle5[stb][en1b][enb];
      EXPAND1(energy, (st + 1, en - 2, DP_U2));
      // (.<   ><   >.) Terminal mismatch
      energy = base_and_branch + arr[st + 2][en - 2][DP_U2] + g_terminal[stb][st1b][en1b][enb];
      EXPAND1(energy, (st + 2, en - 2, DP_U2));

      for (int piv = st + constants::HAIRPIN_MIN_SZ + 2; piv < en - constants::HAIRPIN_MIN_SZ - 2; ++piv) {
        base_t pl1b = r[piv - 1], plb = r[piv], prb = r[piv + 1], pr1b = r[piv + 2];

        // (.(   )   .) Left outer coax - P
        auto outer_coax = energy::MismatchCoaxial(stb, st1b, en1b, enb);
        energy = base_and_branch + arr[st + 2][piv][DP_P] + g_multiloop_hack_b +
            energy::AuGuPenalty(st2b, plb) + arr[piv + 1][en - 2][DP_U] + outer_coax;
        EXPAND2(energy, (st + 2, piv, DP_P), (piv + 1, en - 2, DP_U));

        // (.   (   ).) Right outer coax
        energy = base_and_branch + arr[st + 2][piv][DP_U] + g_multiloop_hack_b +
            energy::AuGuPenalty(prb, en2b) + arr[piv + 1][en - 2][DP_P] + outer_coax;
        EXPAND2(energy, (st + 2, piv, DP_U), (piv + 1, en - 2, DP_P));

        // (.(   ).   ) Left right coax
        energy = base_and_branch + arr[st + 2][piv - 1][DP_P] + g_multiloop_hack_b +
            energy::AuGuPenalty(st2b, pl1b) + arr[piv + 1][en - 1][DP_U] +
            energy::MismatchCoaxial(pl1b, plb, st1b, st2b);
        EXPAND2(energy, (st + 2, piv - 1, DP_P), (piv + 1, en - 1, DP_U));

        // (   .(   ).) Right left coax
        energy = base_and_branch + arr[st + 1][piv][DP_U] + g_multiloop_hack_b +
            energy::AuGuPenalty(pr1b, en2b) + arr[piv + 2][en - 2][DP_P] +
            energy::MismatchCoaxial(en2b, en1b, prb, pr1b);
        EXPAND2(energy, (st + 1, piv, DP_U), (piv + 2, en - 2, DP_P));

        // ((   )   ) Left flush coax
        energy = base_and_branch + arr[st + 1][piv][DP_P] +
            g_multiloop_hack_b + energy::AuGuPenalty(st1b, plb) +
            arr[piv + 1][en - 1][DP_U] + g_stack[stb][st1b][plb][enb];
        EXPAND2(energy, (st + 1, piv, DP_P), (piv + 1, en - 1, DP_U));

        // (   (   )) Right flush coax
        energy = base_and_branch + arr[st + 1][piv][DP_U] +
            g_multiloop_hack_b + energy::AuGuPenalty(prb, en1b) +
            arr[piv + 1][en - 1][DP_P] + g_stack[stb][prb][en1b][enb];
        EXPAND2(energy, (st + 1, piv, DP_U), (piv + 1, en - 1, DP_P));
      }
    } else {
      // Left unpaired. Either DP_U or DP_U2.
      if (st + 1 < en && (a == DP_U || a == DP_U2)) {
        energy = base_energy + arr[st + 1][en][a];
        EXPAND1(energy, (st + 1, en, a));
      }

      // Pair here.
      for (int piv = st + constants::HAIRPIN_MIN_SZ + 1; piv <= en; ++piv) {
        //   (   .   )<   (
        // stb pl1b pb   pr1b
        auto pb = r[piv], pl1b = r[piv - 1];
        // baseAB indicates A bases left unpaired on the left, B bases left unpaired on the right.
        auto base00 = arr[st][piv][DP_P] + energy::AuGuPenalty(stb, pb) + g_multiloop_hack_b;
        auto base01 = arr[st][piv - 1][DP_P] + energy::AuGuPenalty(stb, pl1b) + g_multiloop_hack_b;
        auto base10 = arr[st + 1][piv][DP_P] + energy::AuGuPenalty(st1b, pb) + g_multiloop_hack_b;
        auto base11 = arr[st + 1][piv - 1][DP_P] + energy::AuGuPenalty(st1b, pl1b) + g_multiloop_hack_b;

        // Min is for either placing another unpaired or leaving it as nothing.
        // If we're at U2, don't allow leaving as nothing.
        auto right_unpaired = arr[piv + 1][en][DP_U];
        if (a != DP_U2)
          right_unpaired = std::min(right_unpaired, 0);

        // Check a == U_RCOAX:
        // (   ).<( ** ). > Right coax backward
        if (a == DP_U_RCOAX) {
          if (st > 0) {
            energy = base_energy + base01 + energy::MismatchCoaxial(pl1b, pb, r[st - 1], stb);
            EXPAND1(energy, (st, piv - 1, DP_P));
            EXPAND2(energy + arr[piv + 1][en][DP_U], (st, piv - 1, DP_P), (piv + 1, en, DP_U));
          }
          continue;
        }
        // From here on, a must be U, U2, U_WC, or U_GU.

        // (   )<   > - U, U2, U_WC?, U_GU?
        // Make sure we don't form any branches that are not the right type of pair.
        if ((a != DP_U_WC || IsWatsonCrick(stb, pb)) && (a != DP_U_GU || IsGu(stb, pb))) {
          energy = base_energy + base00;
          // Case where the right doesn't have a branch - for U, U_WC, U_GU (not for U2).
          if (a != DP_U2)
            EXPAND1(energy, (st, piv, DP_P));
          // Case where the right has have a branch - for U, U2, U_WC, U_GU.
          EXPAND2(energy + arr[piv + 1][en][DP_U], (st, piv, DP_P), (piv + 1, en, DP_U));
        }

        // The rest of the cases are for U and U2.
        if (a != DP_U && a != DP_U2)
          continue;

        // (   )3<   > 3' - U, U2
        energy = base_energy + base01 + g_dangle3[pl1b][pb][stb];
        EXPAND1(energy, (st, piv - 1, DP_P));
        EXPAND2(energy + arr[piv + 1][en][DP_U], (st, piv - 1, DP_P), (piv + 1, en, DP_U));

        // 5(   )<   > 5' - U, U2
        energy = base_energy + base10 + g_dangle5[pb][stb][st1b];
        EXPAND1(energy, (st + 1, piv, DP_P));
        EXPAND2(energy + arr[piv + 1][en][DP_U], (st + 1, piv, DP_P), (piv + 1, en, DP_U));

        // .(   ).<   > Terminal mismatch - U, U2
        energy = base_energy + base11 + g_terminal[pl1b][pb][stb][st1b];
        EXPAND1(energy, (st + 1, piv - 1, DP_P));
        EXPAND2(energy + arr[piv + 1][en][DP_U], (st + 1, piv - 1, DP_P), (piv + 1, en, DP_U));

        // .(   ).<(   ) > Left coax - U, U2
        energy = base_energy + base11 + energy::MismatchCoaxial(pl1b, pb, stb, st1b);
        EXPAND2(energy + arr[piv + 1][en][DP_U_WC], (st + 1, piv - 1, DP_P), (piv + 1, en, DP_U_WC));
        EXPAND2(energy + arr[piv + 1][en][DP_U_GU], (st + 1, piv - 1, DP_P), (piv + 1, en, DP_U_GU));

        // (   ).<(   ). > Right coax forward - U, U2
        energy = base_energy + base01 + arr[piv + 1][en][DP_U_RCOAX];
        EXPAND2(energy, (st, piv - 1, DP_P), (piv + 1, en, DP_U_RCOAX));

        // There has to be remaining bases to even have a chance at these cases.
        if (piv < en) {
          auto pr1b = r[piv + 1];
          // (   )<(   ) > Flush coax - U, U2
          energy = base_energy + base00 + g_stack[pb][pr1b][pr1b ^ 3][stb] + arr[piv + 1][en][DP_U_WC];
          EXPAND2(energy, (st, piv, DP_P), (piv + 1, en, DP_U_WC));

          if (pr1b == G || pr1b == U) {
            energy = base_energy + base00 + g_stack[pb][pr1b][pr1b ^ 1][stb] + arr[piv + 1][en][DP_U_GU];
            EXPAND2(energy, (st, piv, DP_P), (piv + 1, en, DP_U_GU));
          }
        }
      }
    }
  }

  std::vector<computed_t> ret;
  for (const auto& struc : finished) {
    for (const auto& meme : struc.not_yet_expanded) {
      printf("%d %d %d\n", std::get<0>(meme), std::get<1>(meme), std::get<2>(meme));
    }
    assert(struc.not_yet_expanded.empty());
    std::vector<std::pair<int, int>> base_pairs;
    for (const auto& hist : struc.history) {
      if (std::get<1>(hist) >= 0 && std::get<2>(hist) == DP_P)
        base_pairs.emplace_back(std::get<0>(hist), std::get<1>(hist));
    }
    // TODO add ctds
    //ret.push_back({r, parsing::BasePairListToPairs(base_pairs, r.size()), struc.energy});
  }
  return ret;
}

#undef PRUNE_INSERT
#undef EXPAND0
#undef EXPAND1
#undef EXPAND2
#undef UNPACK

}
}
