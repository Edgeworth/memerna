// Copyright 2016 Eliot Courtney.
#include "compute/energy/energy.h"

#include <algorithm>
#include <cassert>

#include "compute/energy/internal.h"
#include "compute/energy/structure.h"
#include "model/ctd.h"

namespace mrna::energy {

using internal::AddBranchCtdsToComputed;
using internal::BranchCtd;
using internal::ComputeOptimalCtd;
using internal::GetBranchCtdsFromComputed;

Energy MultiloopEnergy(Computed& computed, bool compute_ctds, int st, int en,
    std::deque<int>& branches, const EnergyModel& em, std::unique_ptr<Structure>* sstruc) {
  const auto& r = computed.r;
  const auto& s = computed.s;
  const bool exterior_loop = s[st] != en;
  Energy energy = 0;

  std::unique_ptr<MultiLoopStructure> struc = nullptr;
  if (sstruc) {
    struc = std::make_unique<MultiLoopStructure>(st, en);
    if (exterior_loop) struc->AddNote("exterior loop");
  }

  // Add AUGU penalties.
  int num_unpaired = 0;
  for (auto branch_st : branches) {
    num_unpaired += s[branch_st] - branch_st + 1;

    if (IsAuGu(r[branch_st], r[s[branch_st]])) {
      if (struc)
        struc->AddNote(
            "%de - opening AU/GU penalty at %d %d", em.augu_penalty, branch_st, s[branch_st]);
      energy += em.augu_penalty;
    }
  }
  num_unpaired = en - st - 1 - num_unpaired + static_cast<int>(exterior_loop) * 2;
  if (struc) struc->AddNote("Unpaired: %d, Branches: %zu", num_unpaired, branches.size() + 1);

  BranchCtd branch_ctds;
  Energy ctd_energy = 0;
  if (exterior_loop) {
    // No initiation for the exterior loop.
    if (compute_ctds) {
      ctd_energy = ComputeOptimalCtd(computed.r, computed.s, em, branches, true, branch_ctds);
      AddBranchCtdsToComputed(computed, branches, branch_ctds);
    } else {
      ctd_energy = GetBranchCtdsFromComputed(computed, em, branches, branch_ctds);
    }
  } else {
    if (IsAuGu(r[st], r[en])) {
      if (struc) struc->AddNote("%de - closing AU/GU penalty at %d %d", em.augu_penalty, st, en);
      energy += em.augu_penalty;
    }
    Energy initiation = em.MultiloopInitiation(static_cast<int>(branches.size() + 1));
    if (struc) struc->AddNote("%de - initiation", initiation);
    energy += initiation;

    if (compute_ctds) {
      BranchCtd config_ctds[4] = {};
      std::pair<Energy, int> config_energies[4] = {};
      branches.push_front(en);
      config_energies[0] = {
          ComputeOptimalCtd(computed.r, computed.s, em, branches, true, config_ctds[0]), 0};
      config_energies[1] = {
          ComputeOptimalCtd(computed.r, computed.s, em, branches, false, config_ctds[1]), 1};
      branches.pop_front();
      branches.push_back(en);
      config_energies[2] = {
          ComputeOptimalCtd(computed.r, computed.s, em, branches, true, config_ctds[2]), 2};
      // Swap the final branch back to the front because following code expects it.
      config_ctds[2].push_front(config_ctds[2].back());
      config_ctds[2].pop_back();
      config_energies[3] = {
          ComputeOptimalCtd(computed.r, computed.s, em, branches, false, config_ctds[3]), 3};
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

  if (struc) {
    struc->AddNote("%de - ctd", ctd_energy);
    if (!exterior_loop) {
      struc->AddNote("%de - outer loop stacking - %s", branch_ctds[0].second,
          energy::CtdToName(branch_ctds[0].first));
      branch_ctds.pop_front();
    }
    for (const auto& ctd : branch_ctds) struc->AddCtd(ctd.first, ctd.second);
    // Give the pointer back.
    *sstruc = std::move(struc);
  }

  return energy;
}

Energy ComputeSubstructureEnergy(Computed& computed, bool compute_ctds, int st, int en,
    const EnergyModel& em, std::unique_ptr<Structure>* struc) {
  const auto& r = computed.r;
  const auto& s = computed.s;
  const bool exterior_loop = s[st] != en;
  assert(en >= st);
  Energy energy = 0;

  // Look for branches inside.
  std::deque<int> branches;
  for (int i = st; i <= en; ++i) {
    int pair = s[i];
    assert(pair <= en && (pair == -1 || s[pair] == i));
    if (!(i == st && pair == en) && !(i == en && pair == st) && pair != -1) {
      branches.push_back(i);
      // Skip ahead.
      i = pair;
    }
  }

  if (exterior_loop || branches.size() >= 2) {
    // Multiloop.
    energy += MultiloopEnergy(computed, compute_ctds, st, en, branches, em, struc);
  } else if (branches.empty()) {
    // Hairpin loop.
    assert(en - st - 1 >= 3);
    energy += em.Hairpin(r, st, en, struc);
  } else if (branches.size() == 1) {
    const int loop_st = branches.front(), loop_en = s[branches.front()];
    energy += em.TwoLoop(r, st, en, loop_st, loop_en, struc);
  }

  if (struc) (*struc)->set_self_energy(energy);
  // Add energy from children.
  for (auto i : branches) {
    if (struc) {
      std::unique_ptr<Structure> structure;
      energy += ComputeSubstructureEnergy(computed, compute_ctds, i, s[i], em, &structure);
      (*struc)->AddBranch(std::move(structure));
    } else {
      energy += ComputeSubstructureEnergy(computed, compute_ctds, i, s[i], em, nullptr);
    }
  }
  if (struc) (*struc)->set_total_energy(energy);

  return energy;
}

Computed ComputeEnergy(const Primary& r, const Secondary& s, const EnergyModel& em,
    std::unique_ptr<Structure>* struc) {
  Computed computed(r, s);
  return ComputeEnergyWithCtds(computed, em, true, struc);
}

Computed ComputeEnergyWithCtds(const Computed& computed, const EnergyModel& em, bool compute_ctds,
    std::unique_ptr<Structure>* struc) {
  auto computed_copy = computed;
  const auto& r = computed_copy.r;
  const auto& s = computed_copy.s;
  Energy energy = ComputeSubstructureEnergy(
      computed_copy, compute_ctds, 0, static_cast<int>(r.size()) - 1, em, struc);
  if (s[0] == static_cast<int>(r.size() - 1) && IsAuGu(r[0], r[s[0]])) {
    energy += em.augu_penalty;
    if (struc) {
      (*struc)->AddNote("%de - top level AU/GU penalty", em.augu_penalty);
      (*struc)->set_self_energy((*struc)->self_energy() + em.augu_penalty);  // Gross.
      (*struc)->set_total_energy((*struc)->total_energy() + em.augu_penalty);  // Gross.
    }
  }
  computed_copy.energy = energy;
  return computed_copy;
}

}  // namespace mrna::energy
