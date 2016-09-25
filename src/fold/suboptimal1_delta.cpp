#include <stack>
#include "fold/suboptimal1_base.h"
#include "fold/suboptimal1_delta.h"

namespace memerna {
namespace fold {
namespace internal {

void Suboptimal1Delta::Run(std::function<void(const computed_t&)> fn) {
  // TODO merge max-delta functionality into this, with repeated dfs? at least use for sorted output
  memset(gp.data(), -1, gp.size());
  memset(gctd.data(), CTD_NA, gctd.size());

  // General idea is perform a dfs of the expand tree. Keep track of the current partial structures
  // and energy. Also keep track of what is yet to be expanded. Each node is either a terminal,
  // or leads to one expansion (either from unexpanded, or from expanding itself) - if there is
  // another expansion it is put on the unexpanded list. Everything in unexpanded never affects
  // the CTDs or energy of the current state - that is rolled into the modification when that
  // unexpanded is originally generated.

  energy_t energy = 0;
  q.push({0, {0, -1, EXT}, false});
  while (q.size()) {
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

    // If we ran out of expansions, or the next expansion would take us over the delta limit
    // we are done with this node.
    if (s.idx == int(exps.size()) || exps[s.idx].energy + energy > delta) {
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
        computed_t tmp_computed = {{
            std::move(gr), std::move(gp)}, std::move(gctd), energy + gext[0][EXT]};
        fn(tmp_computed);
        // Move everything back
        gr = std::move(tmp_computed.s.r);
        gp = std::move(tmp_computed.s.p);
        gctd = std::move(tmp_computed.base_ctds);
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
}

}
}
}
