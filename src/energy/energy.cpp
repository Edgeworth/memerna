#include <cstdio>
#include <cmath>
#include <memory>
#include <algorithm>
#include "energy/energy.h"
#include "energy/structure.h"

namespace memerna {
namespace energy {

using namespace internal;

energy_t MultiloopEnergy(computed_t& computed, bool compute_ctds,
    int st, int en, std::deque<int>& branches, const EnergyModel& em,
    std::unique_ptr<Structure>* ss) {
  const auto& r = computed.s.r;
  const auto& p = computed.s.p;
  const bool exterior_loop = st == 0 && en == int(r.size() - 1) && p[st] != en;
  energy_t energy = 0;

  std::unique_ptr<MultiLoopStructure> s = nullptr;
  if (ss) {
    s = std::make_unique<MultiLoopStructure>(st, en);
    if (exterior_loop) s->AddNote("exterior loop");
  }

  // Add AUGU penalties.
  int num_unpaired = 0;
  for (auto branch_st : branches) {
    num_unpaired += p[branch_st] - branch_st + 1;

    if (IsAuGu(r[branch_st], r[p[branch_st]])) {
      if (s) s->AddNote("%de - opening AU/GU penalty at %d %d", em.augu_penalty, branch_st, p[branch_st]);
      energy += em.augu_penalty;
    }
  }
  num_unpaired = en - st - 1 - num_unpaired + exterior_loop * 2;
  if (s) s->AddNote("Unpaired: %d, Branches: %zu", num_unpaired, branches.size() + 1);

  branch_ctd_t branch_ctds;
  energy_t ctd_energy = 0;
  if (exterior_loop) {
    // No initiation for the exterior loop.
    if (compute_ctds) {
      ctd_energy = ComputeOptimalCtd(computed.s, em, branches, true, branch_ctds);
      AddBranchCtdsToComputed(computed, branches, branch_ctds);
    } else {
      ctd_energy = GetBranchCtdsFromComputed(computed, em, branches, branch_ctds);
    }
  } else {
    if (IsAuGu(r[st], r[en])) {
      if (s) s->AddNote("%de - closing AU/GU penalty at %d %d", em.augu_penalty, st, en);
      energy += em.augu_penalty;
    }
    energy_t initiation = em.MultiloopInitiation(int(branches.size() + 1));
    if (s) s->AddNote("%de - initiation", initiation);
    energy += initiation;

    if (compute_ctds) {
      branch_ctd_t config_ctds[4] = {};
      std::pair<energy_t, int> config_energies[4] = {};
      branches.push_front(en);
      config_energies[0] = {ComputeOptimalCtd(computed.s, em, branches, true, config_ctds[0]), 0};
      config_energies[1] = {ComputeOptimalCtd(computed.s, em, branches, false, config_ctds[1]), 1};
      branches.pop_front();
      branches.push_back(en);
      config_energies[2] = {ComputeOptimalCtd(computed.s, em, branches, true, config_ctds[2]), 2};
      // Swap the final branch back to the front because following code expects it.
      config_ctds[2].push_front(config_ctds[2].back());
      config_ctds[2].pop_back();
      config_energies[3] = {ComputeOptimalCtd(computed.s, em, branches, false, config_ctds[3]), 3};
      config_ctds[3].push_front(config_ctds[3].back());
      config_ctds[3].pop_back();
      branches.pop_back();
      std::sort(config_energies, config_energies + 4);
      branch_ctds = config_ctds[config_energies[0].second];
      ctd_energy = config_energies[0].first;

      // Write the optimal ctds to computed.
      branches.push_front(en);
      AddBranchCtdsToComputed(computed, branches, branch_ctds);
      branches.pop_front();
    } else {
      branches.push_front(en);
      ctd_energy = GetBranchCtdsFromComputed(computed, em, branches, branch_ctds);
      branches.pop_front();
    }
  }
  energy += ctd_energy;

  if (s) {
    s->AddNote("%de - ctd", ctd_energy);
    if (!exterior_loop) {
      s->AddNote("%de - outer loop stacking - %s", branch_ctds[0].second,
          energy::CtdToName(branch_ctds[0].first));
      branch_ctds.pop_front();
    }
    for (const auto& ctd : branch_ctds)
      s->AddCtd(ctd.first, ctd.second);
    // Give the pointer back.
    ss->reset(s.release());
  }

  return energy;
}

energy_t ComputeEnergyInternal(computed_t& computed, bool compute_ctds,
    int st, int en, const EnergyModel& em, std::unique_ptr<Structure>* s) {
  const auto& r = computed.s.r;
  const auto& p = computed.s.p;
  assert(en >= st);
  energy_t energy = 0;

  // Look for branches inside.
  std::deque<int> branches;
  for (int i = st; i <= en; ++i) {
    int pair = p[i];
    assert(pair <= en);
    if (!(i == st && pair == en) && !(i == en && pair == st) && pair != -1) {
      branches.push_back(i);
      // Skip ahead.
      i = pair;
    }
  }

  // We're in the exterior loop if we were called with the entire RNA and there's no match on the very ends that takes
  // us out of the exterior loop.
  const bool exterior_loop = st == 0 && en == int(r.size() - 1) && p[st] != en;
  if (exterior_loop || branches.size() >= 2) {
    // Multiloop.
    energy += MultiloopEnergy(computed, compute_ctds, st, en, branches, em, s);
  } else if (branches.size() == 0) {
    // Hairpin loop.
    assert(en - st - 1 >= 3);
    energy += em.Hairpin(r, st, en, s);
  } else if (branches.size() == 1) {
    const int loop_st = branches.front(), loop_en = p[branches.front()];
    energy += em.TwoLoop(r, st, en, loop_st, loop_en, s);
  }

  if (s) (*s)->SetSelfEnergy(energy);
  // Add energy from children.
  for (auto i : branches) {
    if (s) {
      std::unique_ptr<Structure> structure;
      energy += ComputeEnergyInternal(computed, compute_ctds, i, p[i], em, &structure);
      (*s)->AddBranch(std::move(structure));
    } else {
      energy += ComputeEnergyInternal(computed, compute_ctds, i, p[i], em, nullptr);
    }
  }
  if (s) (*s)->SetTotalEnergy(energy);

  return energy;
}


computed_t ComputeEnergy(const secondary_t& secondary,
    const EnergyModel& em, std::unique_ptr<Structure>* s) {
  computed_t computed(secondary);
  return ComputeEnergyWithCtds(computed, em, true, s);
}

computed_t ComputeEnergyWithCtds(const computed_t& computed,
    const EnergyModel& em, bool compute_ctds, std::unique_ptr<Structure>* s) {
  auto computed_copy = computed;
  const auto& r = computed_copy.s.r;
  const auto& p = computed_copy.s.p;
  energy_t energy = ComputeEnergyInternal(computed_copy, compute_ctds, 0, (int) r.size() - 1, em, s);
  if (p[0] == int(r.size() - 1) && IsAuGu(r[0], r[p[0]])) {
    energy += em.augu_penalty;
    if (s) {
      (*s)->AddNote("%de - top level AU/GU penalty", em.augu_penalty);
      (*s)->SetSelfEnergy((*s)->GetSelfEnergy() + em.augu_penalty);  // Gross.
      (*s)->SetTotalEnergy((*s)->GetTotalEnergy() + em.augu_penalty);  // Gross.
    }
  }
  computed_copy.energy = energy;
  return computed_copy;
}


}
}
