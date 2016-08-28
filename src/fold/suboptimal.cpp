#include <set>
#include "suboptimal.h"
#include "parsing.h"

namespace memerna {
namespace fold {

struct suboptimal_node_t {

  // Since nodes form a tree, can save memory on these two vectors.
  // Will also save time in the comparison.
  std::vector<std::tuple<int, int, int>> not_yet_expanded;
  std::vector<std::pair<int, int>> base_pairs;
  energy_t energy;  // Stores the minimum energy this state could have.

  bool operator<(const suboptimal_node_t& o) const {
    if (energy != o.energy) return energy < o.energy;
    if (not_yet_expanded != o.not_yet_expanded) return not_yet_expanded < o.not_yet_expanded;
    assert(false);  // This shouldn't happen.
    if (base_pairs != o.base_pairs) return base_pairs < o.base_pairs;
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
std::vector<folded_rna_t> SuboptimalTraceback(
    energy_t max_energy, int max_structures,
    const array3d_t<energy_t, DP_SIZE>& arr,
    const array2d_t<energy_t, EXT_SIZE>& exterior) {
  // Basic idea of suboptimal traceback is look at all possible choices from a state, and expand just one of them.
  // Fully expanding one of them means there will be no duplicates in the tree.
  // Cull the ones not inside the window or when we have more than |max_structures|.
  verify_expr(max_structures > 0, "must request at least one structure");
  int N = int(r.size());
  std::set<suboptimal_node_t> finished;
  std::set<suboptimal_node_t> q;

  suboptimal_node_t newnode = {{std::make_tuple(0, -1, EXT)}, {}, exterior[0][EXT]};
  q.insert(newnode);
  while (!q.empty()) {
    printf("meme: %zu\n", q.size());
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
    int st, en, a;
    std::tie(st, en, a) = to_expand;
    printf("meme: %zu %d %d %d\n", q.size(), st, en, a);

    // Initialise, since we only make small modifications to it.
    newnode = node;

    // Exterior loop
    if (en == -1) {
      // We try replacing what we do at (st, a) with a bunch of different cases, so we use this energy as a base.
      energy_t base_energy = node.energy - exterior[st][a];
      // Base case: do nothing.
      newnode.energy = base_energy;
      PRUNE_INSERT(q, newnode);
      // Case: No pair starting here (for EXT only)
      if (a == EXT) {
        newnode.not_yet_expanded.emplace_back(st + 1, -1, EXT);
        newnode.energy = base_energy + exterior[st + 1][EXT];
        PRUNE_INSERT(q, newnode);
        newnode.not_yet_expanded.pop_back();
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
        // Most take this.
        newnode.not_yet_expanded.emplace_back(en + 1, -1, EXT);

        // (   ).<( * ). > Right coax backward
        if (st > 0 && a == EXT_RCOAX) {
          newnode.not_yet_expanded.emplace_back(st, en - 1, DP_P);
          newnode.energy = base_energy + base01 + energy::MismatchCoaxial(
              en1b, enb, r[st - 1], stb) + exterior[en + 1][EXT];
          PRUNE_INSERT(q, newnode);
          newnode.not_yet_expanded.pop_back();
        }

        // Cases for EXT, EXT_WC, EXT_GU.
        if (a == EXT || (a == EXT_WC && IsWatsonCrick(stb, enb)) || (a == EXT_GU && IsGu(stb, enb))) {
          // (   )<   >
          newnode.not_yet_expanded.emplace_back(st, en, DP_P);
          newnode.energy = base_energy + base00 + exterior[en + 1][EXT];
          PRUNE_INSERT(q, newnode);
          newnode.not_yet_expanded.pop_back();

          // (   )3<   > 3'
          newnode.not_yet_expanded.emplace_back(st, en - 1, DP_P);
          newnode.energy = base_energy + base01 + g_dangle3_e[en1b][enb][stb] + exterior[en + 1][EXT];
          PRUNE_INSERT(q, newnode);
          newnode.not_yet_expanded.pop_back();
        }

        // Everything after this is only for EXT.
        if (a != EXT) continue;

        // 5(   )<   > 5'
        newnode.not_yet_expanded.emplace_back(st + 1, en, DP_P);
        newnode.energy = base_energy + base10 + g_dangle5_e[enb][stb][st1b] + exterior[en + 1][EXT];
        PRUNE_INSERT(q, newnode);
        newnode.not_yet_expanded.pop_back();

        // .(   ).<   > Terminal mismatch
        newnode.not_yet_expanded.emplace_back(st + 1, en - 1, DP_P);
        newnode.energy = base_energy + base11 + g_terminal[en1b][enb][stb][st1b] + exterior[en + 1][EXT];
        PRUNE_INSERT(q, newnode);
        newnode.not_yet_expanded.pop_back();

        // .(   ).<(   ) > Left coax
        newnode.not_yet_expanded.pop_back();  // Remove EXT
        newnode.not_yet_expanded.emplace_back(st + 1, en - 1, DP_P);
        auto val = base_energy + base11 + energy::MismatchCoaxial(en1b, enb, stb, st1b);

        newnode.not_yet_expanded.emplace_back(en + 1, -1, EXT_GU);
        newnode.energy = val + exterior[en + 1][EXT_GU];
        PRUNE_INSERT(q, newnode);
        newnode.not_yet_expanded.pop_back();

        newnode.not_yet_expanded.emplace_back(en + 1, -1, EXT_WC);
        newnode.energy = val + exterior[en + 1][EXT_WC];
        PRUNE_INSERT(q, newnode);
        newnode.not_yet_expanded.pop_back();
        newnode.not_yet_expanded.pop_back();  // Also remove EXT_WC.

        // (   ).<(   ). > Right coax forward
        newnode.not_yet_expanded.emplace_back(en + 1, -1, EXT_RCOAX);
        newnode.not_yet_expanded.emplace_back(st, en - 1, DP_P);
        newnode.energy = base_energy + base01 + exterior[en + 1][EXT_RCOAX];
        PRUNE_INSERT(q, newnode);
        newnode.not_yet_expanded.pop_back();

        if (en < N - 1) {
          // (   )<(   ) > Flush coax
          auto enrb = r[en + 1];
          newnode.not_yet_expanded.emplace_back(st, en, DP_P);
          newnode.not_yet_expanded.emplace_back(en + 1, -1, EXT_WC);
          newnode.energy = base_energy + base00 + g_stack[enb][enrb][enrb ^ 3][stb] + exterior[en + 1][EXT_WC];
          PRUNE_INSERT(q, newnode);
          newnode.not_yet_expanded.pop_back();

          if (enrb == G || enrb == U) {
            newnode.not_yet_expanded.emplace_back(en + 1, -1, EXT_GU);
            newnode.energy = base_energy + base00 + g_stack[enb][enrb][enrb ^ 1][stb] + exterior[en + 1][EXT_GU];
            PRUNE_INSERT(q, newnode);
            newnode.not_yet_expanded.pop_back();
          }
        }
      }
      // Finished exterior loop, don't do anymore.
      continue;
    }

    // Normal stuff
    if (a == DP_P) {
      newnode.base_pairs.emplace_back(st, en);
      PRUNE_INSERT(q, newnode);
    }
  }

  std::vector<folded_rna_t> ret;
  for (const auto& struc : finished) {
    for (const auto& meme : struc.not_yet_expanded) {
      printf("%d %d %d\n", std::get<0>(meme), std::get<1>(meme), std::get<2>(meme));
    }
    assert(struc.not_yet_expanded.empty());
    ret.push_back({r, parsing::BasePairListToPairs(struc.base_pairs, r.size()), struc.energy});
  }
  return ret;
}
#undef PRUNE_INSERT
}
}
