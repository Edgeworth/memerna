// Copyright 2016 Eliot Courtney.
#include "compute/energy/t04/model.h"

#include <fmt/core.h>

#include <cassert>
#include <deque>
#include <memory>
#include <utility>

#include "compute/energy/energy.h"
#include "compute/energy/structure.h"
#include "model/base.h"
#include "model/ctd.h"
#include "model/primary.h"

namespace mrna::erg::t04 {

Energy Model::SubEnergyInternal(const Primary& r, const Secondary& s, int st, int en,
    bool use_given_ctds, Ctds* ctd, std::unique_ptr<Structure>* struc) const {
  assert(en >= st);
  const bool exterior_loop = s[st] != en;
  Energy energy = ZERO_E;

  // Look for branches inside.
  std::deque<int> branches;
  for (int i = st; i <= en; ++i) {
    const int pair = s[i];
    assert(pair <= en && (pair == -1 || s[pair] == i));
    if ((i != st || pair != en) && (i != en || pair != st) && pair != -1) {
      branches.push_back(i);
      // Skip ahead.
      i = pair;
    }
  }

  if (exterior_loop || branches.size() >= 2) {
    // Multiloop.
    energy += MultiloopEnergy(r, s, st, en, &branches, use_given_ctds, ctd, struc);
  } else if (branches.empty()) {
    // Hairpin loop.
    assert(en - st - 1 >= 3);
    energy += Hairpin(r, st, en, struc);
  } else {
    assert(branches.size() == 1);
    const int loop_st = branches.front();
    const int loop_en = s[branches.front()];
    energy += TwoLoop(r, st, en, loop_st, loop_en, struc);
  }

  if (struc) (*struc)->set_self_energy(energy);
  // Add energy from children.
  for (auto i : branches) {
    if (struc) {
      std::unique_ptr<Structure> sstruc;
      energy += SubEnergyInternal(r, s, i, s[i], use_given_ctds, ctd, &sstruc);
      (*struc)->AddBranch(std::move(sstruc));
    } else {
      energy += SubEnergyInternal(r, s, i, s[i], use_given_ctds, ctd, nullptr);
    }
  }
  if (struc) (*struc)->set_total_energy(energy);

  return energy;
}

// If (st, en) is not paired, treated as an exterior loop.
// If |ctd| is non-null, use the given ctds.
EnergyResult Model::SubEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd, int st,
    int en, bool build_structure) const {
  const bool use_given_ctds = given_ctd;
  auto ctd = use_given_ctds ? Ctds(*given_ctd) : Ctds(r.size());

  std::unique_ptr<Structure> struc;
  auto energy =
      SubEnergyInternal(r, s, st, en, use_given_ctds, &ctd, build_structure ? &struc : nullptr);
  return {energy, std::move(ctd), std::move(struc)};
}

EnergyResult Model::TotalEnergy(
    const Primary& r, const Secondary& s, const Ctds* given_ctd, bool build_structure) const {
  auto res = SubEnergy(r, s, given_ctd, 0, static_cast<int>(r.size()) - 1, build_structure);
  if (s[0] == static_cast<int>(r.size() - 1) && IsAuPair(r[0], r[s[0]])) {
    res.energy += au_penalty;
    if (res.struc) {
      res.struc->AddNote("{}e - top level AU penalty", au_penalty);
      res.struc->set_self_energy(res.struc->self_energy() + au_penalty);
      res.struc->set_total_energy(res.struc->total_energy() + au_penalty);
    }
  }
  if (s[0] == static_cast<int>(r.size() - 1) && IsGuPair(r[0], r[s[0]])) {
    res.energy += gu_penalty;
    if (res.struc) {
      res.struc->AddNote("{}e - top level GU penalty", gu_penalty);
      res.struc->set_self_energy(res.struc->self_energy() + gu_penalty);
      res.struc->set_total_energy(res.struc->total_energy() + gu_penalty);
    }
  }
  return res;
}

}  // namespace mrna::erg::t04
