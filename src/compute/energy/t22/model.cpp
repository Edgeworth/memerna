// Copyright 2022 Eliot Courtney.
#include "compute/energy/t22/model.h"

#include <fmt/core.h>

#include <cassert>
#include <deque>
#include <memory>
#include <random>
#include <string>
#include <utility>

#include "compute/energy/common/model.h"
#include "compute/energy/common/parse.h"
#include "compute/energy/common/t04like/model_mixin.h"
#include "compute/energy/energy.h"
#include "compute/energy/structure.h"
#include "model/base.h"
#include "model/ctd.h"
#include "model/primary.h"

namespace mrna::erg::t22 {

Energy Model::SubEnergyInternal(const Primary& r, const Secondary& s, int st, int en, int stack_st,
    int stack_en, bool use_given_ctds, Ctds* ctd, std::unique_ptr<Structure>* struc) const {
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

  bool continuous = false;
  if (exterior_loop || branches.size() >= 2) {
    // Multiloop.
    energy += MultiloopEnergy(r, s, st, en, &branches, use_given_ctds, ctd, struc);

    // Current stack is terminated.
    if (stack_st != -1 && stack_st != st)
      energy += StackPenalty(r, s, stack_st, stack_en, st, en, struc);
    stack_st = -1;
    stack_en = -1;
  } else if (branches.empty()) {
    // Hairpin loop.
    assert(en - st - 1 >= 3);
    energy += Hairpin(r, st, en, struc);

    // Current stack is terminated.
    if (stack_st != -1 && stack_st != st)
      energy += StackPenalty(r, s, stack_st, stack_en, st, en, struc);
    stack_st = -1;
    stack_en = -1;
  } else {
    assert(branches.size() == 1);
    const int loop_st = branches.front();
    const int loop_en = s[branches.front()];
    energy += TwoLoop(r, st, en, loop_st, loop_en, struc);

    // Current stack is terminated if it's not continuous, otherwise it's extended.
    if (!IsContinuous(st, en, loop_st, loop_en)) {
      if (stack_st != -1 && stack_st != st)
        energy += StackPenalty(r, s, stack_st, stack_en, st, en, struc);
      stack_st = -1;
      stack_en = -1;
    } else {
      continuous = true;
    }
  }

  if (struc) (*struc)->set_self_energy(energy);
  // Add energy from children.
  for (auto i : branches) {
    int nst = continuous ? stack_st : i;
    int nen = continuous ? stack_en : s[i];
    if (struc) {
      std::unique_ptr<Structure> sstruc;
      energy += SubEnergyInternal(r, s, i, s[i], nst, nen, use_given_ctds, ctd, &sstruc);
      (*struc)->AddBranch(std::move(sstruc));
    } else {
      energy += SubEnergyInternal(r, s, i, s[i], nst, nen, use_given_ctds, ctd, nullptr);
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

  const bool exterior_loop = s[st] != en;
  int stack_st = exterior_loop ? -1 : st;
  int stack_en = exterior_loop ? -1 : en;

  std::unique_ptr<Structure> struc;
  auto energy = SubEnergyInternal(
      r, s, st, en, stack_st, stack_en, use_given_ctds, &ctd, build_structure ? &struc : nullptr);
  return {energy, std::move(ctd), std::move(struc)};
}

EnergyResult Model::TotalEnergy(
    const Primary& r, const Secondary& s, const Ctds* given_ctd, bool build_structure) const {
  auto res = SubEnergy(r, s, given_ctd, 0, static_cast<int>(r.size()) - 1, build_structure);
  if (s[0] == static_cast<int>(r.size() - 1)) {
    auto extra_energy = ZERO_E;
    if (IsAuPair(r[0], r[s[0]])) {
      extra_energy += au_penalty;
      if (res.struc) res.struc->AddNote("{}e - top level AU penalty", au_penalty);
    }
    if (IsGuPair(r[0], r[s[0]])) {
      extra_energy += gu_penalty;
      if (res.struc) res.struc->AddNote("{}e - top level GU penalty", gu_penalty);
    }
    res.energy += extra_energy;
    if (res.struc) {
      res.struc->set_self_energy(res.struc->self_energy() + extra_energy);
      res.struc->set_total_energy(res.struc->total_energy() + extra_energy);
    }
  }

  return res;
}

void Model::LoadFromDir(const std::string& data_dir) {
  T04ModelMixin::LoadFromDir(data_dir);
  Parse4MapFromFile(data_dir + "/penultimate_stacking.data", penultimate_stack);
}

void Model::LoadRandom(std::mt19937& eng) {
  T04ModelMixin::LoadRandom(eng);

  // penultimate_stack is dependent on the direction, so 180 degree rotations
  // don't have to be the same.
  std::uniform_real_distribution<double> energy_dist(RAND_MIN_ENERGY, RAND_MAX_ENERGY);
  RANDOMISE_DATA(penultimate_stack);
}

}  // namespace mrna::erg::t22
