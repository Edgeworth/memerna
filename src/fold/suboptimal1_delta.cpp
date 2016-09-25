#include <stack>
#include "fold/suboptimal1_base.h"
#include "fold/suboptimal1_delta.h"

namespace memerna {
namespace fold {
namespace internal {

int Suboptimal1Delta::Run(std::function<void(const computed_t&)> fn, bool sorted) {
  // TODO merge max-delta functionality into this, with repeated dfs? at least use for sorted output
  memset(gp.data(), -1, gp.size());
  memset(gctd.data(), CTD_NA, gctd.size());

  // If require sorted output, or limited number of structures (requires sorting).
  if (sorted || max_structures != MAX_STRUCTURES) {
    int num_structures = 0;
    energy_t cur_delta = 0;
    while (num_structures < max_structures && cur_delta != MAX_E && cur_delta <= delta) {
      auto res = RunInternal(fn, cur_delta, true, max_structures - num_structures);
      num_structures += res.first;
      cur_delta = res.second;
    }
    return num_structures;
  }

  return RunInternal(fn, delta, false, MAX_STRUCTURES).first;
}

std::pair<int, int> Suboptimal1Delta::RunInternal(std::function<void(const computed_t&)> fn,
    energy_t cur_delta, bool exact_energy, int structure_limit) {
  // General idea is perform a dfs of the expand tree. Keep track of the current partial structures
  // and energy. Also keep track of what is yet to be expanded. Each node is either a terminal,
  // or leads to one expansion (either from unexpanded, or from expanding itself) - if there is
  // another expansion it is put on the unexpanded list. Everything in unexpanded never affects
  // the CTDs or energy of the current state - that is rolled into the modification when that
  // unexpanded is originally generated.

  int num_structures = 0;
  // Store the smallest energy above cur_delta we see. If we reach our |structure_limit| before
  // finishing, we might not see the smallest one, but it's okay since we won't be called again.
  // Otherwise, we will completely finish, and definitely see it.
  energy_t next_seen = MAX_E;
  energy_t energy = 0;
  q.push({0, {0, -1, EXT}, false});
  while (!q.empty()) {
    auto& s = q.top();
    assert(s.expand.st != -1);

    const auto& exps = GetExpansions(s.expand);
    assert(!exps.empty());  // Must produce at least one expansion: {-1, -1, -1}.

    // Undo previous child's ctds and energy. The pairing is undone by the child.
    // Also remove from unexpanded if the previous child added stuff to it.
    if (s.idx != 0) {
      const auto& pexp = exps[s.idx - 1];
      if (pexp.ctd0.idx != -1) gctd[pexp.ctd0.idx] = CTD_NA;
      if (pexp.ctd1.idx != -1) gctd[pexp.ctd1.idx] = CTD_NA;
      if (pexp.unexpanded.st != -1) unexpanded.pop_back();
      energy -= pexp.energy;
    }

    // Update the next best seen variable
    if (s.idx != int(exps.size()) && exps[s.idx].energy + energy > cur_delta)
      next_seen = std::min(next_seen, exps[s.idx].energy + energy);

    // If we ran out of expansions, or the next expansion would take us over the delta limit
    // we are done with this node.
    if (s.idx == int(exps.size()) || exps[s.idx].energy + energy > cur_delta) {
      // Finished looking at this node, so undo this node's modifications to the global state.
      if (s.expand.en != -1 && s.expand.a == DP_P)
        gp[s.expand.st] = gp[s.expand.en] = -1;
      if (s.should_unexpand)
        unexpanded.push_back(s.expand);
      q.pop();
      continue;  // Done.
    }

    const auto& exp = exps[s.idx++];
    dfs_state_t ns = {0, exp.to_expand, false};
    energy += exp.energy;
    if (exp.to_expand.st == -1) {
      // Can't have an unexpanded without a to_expand. Also can't set ctds or affect energy.
      assert(exp.unexpanded.st == -1);
      assert(exp.ctd0.idx == -1 && exp.ctd1.idx == -1);
      // Use an unexpanded now, if one exists.
      if (unexpanded.empty()) {
        // At a terminal state.
        if (!exact_energy || energy == cur_delta) {
          computed_t tmp_computed = {{
              std::move(gr), std::move(gp)}, std::move(gctd), energy + gext[0][EXT]};
          fn(tmp_computed);
          ++num_structures;
          // Move everything back
          gr = std::move(tmp_computed.s.r);
          gp = std::move(tmp_computed.s.p);
          gctd = std::move(tmp_computed.base_ctds);

          // Hit structure limit.
          if (num_structures == structure_limit)
            return {num_structures, -1};
        }
        continue;  // Done
      } else {
        ns.expand = unexpanded.back();
        unexpanded.pop_back();
        // This node should replace itself into |unexpanded| when its done.
        ns.should_unexpand = true;
      }
    } else {
      // Apply child's modifications to the global state.
      if (exp.ctd0.idx != -1) gctd[exp.ctd0.idx] = exp.ctd0.ctd;
      if (exp.ctd1.idx != -1) gctd[exp.ctd1.idx] = exp.ctd1.ctd;
      if (exp.unexpanded.st != -1)
        unexpanded.push_back(exp.unexpanded);
    }
    if (ns.expand.en != -1 && ns.expand.a == DP_P) {
      gp[ns.expand.st] = ns.expand.en;
      gp[ns.expand.en] = ns.expand.st;
    }
    q.push(ns);
  }
  assert(unexpanded.empty() && energy == 0 &&
      gp == std::vector<int>(gp.size(), -1) &&
      gctd == std::vector<Ctd>(gctd.size(), CTD_NA));
  return {num_structures, next_seen};
}

}
}
}
