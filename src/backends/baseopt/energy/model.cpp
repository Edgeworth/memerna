// Copyright 2016 Eliot Courtney.
#include "backends/baseopt/energy/model.h"

#include <fmt/core.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <deque>
#include <memory>
#include <optional>
#include <utility>

#include "api/energy/energy_cfg.h"
#include "backends/common/base/branch.h"
#include "backends/common/base/model_base.h"
#include "backends/common/branch.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/structure.h"

namespace mrna::md::base::opt {

// Indices are inclusive, include the initiating base pair.
// N.B. This includes an ending AU/GU penalty.
// Rules for hairpin energy:
// 1. Check a lookup table for a hardcoded value.
// 2. Loops with less than 3 bases on the inside are disallowed by base.
// Case 1 -- 3 bases on the inside:
//   Return G_init + penalty if all bases inside are C.
//   The penalty for all C bases for size 3 hairpin loops is special cased.
// Case 2 -- more than 3 bases on the inside:
//   Return G_init plus the following penalties / bonuses:
//   Terminal mismatch energy for st + 1 and en - 1.
//   If the mismatch is UU or GA (not AG), additional bonus
//   If the mismatch is GG, additional bonus.
//   If the pair st, en is GU (not UG), a bonus if st - 1 and st - 2 are both Gs, if they exist.
//   A penalty if all the bases inside are C: A * length + B (A, B specified as part of the energy
//   model).
Energy Model::Hairpin(const Primary& r, int st, int en, std::unique_ptr<Structure>* s) const {
  assert(st < en);
  if (s) *s = std::make_unique<HairpinLoopStructure>(st, en);

  Energy energy = ZERO_E;

  // Apply AU penalty if necessary.
  if (IsAuPair(r[st], r[en])) {
    if (s) (*s)->AddNote("{}e - AU penalty", au_penalty);
    energy += au_penalty;
  }
  if (IsGuPair(r[st], r[en])) {
    if (s) (*s)->AddNote("{}e - GU penalty", gu_penalty);
    energy += gu_penalty;
  }

  std::string seq;
  for (int i = st; i <= en; ++i) seq += BaseToChar(r[i]).value();
  const auto iter = hairpin.find(seq);
  if (iter != hairpin.end()) {
    if (s) (*s)->AddNote("special hairpin");
    return energy + iter->second;
  }

  // Subtract two for the initiating base pair.
  const int length = en - st - 1;
  if (length < 3) return MAX_E;  // Disallowed by base.

  const auto initiation = HairpinInitiation(length);
  energy += initiation;
  if (s) (*s)->AddNote("{}e - initiation", initiation);
  // base says hairpin loops with all C bases inside them are treated specially.
  bool all_c = true;
  for (int i = st + 1; i <= en - 1; ++i) {
    if (r[i] != C) all_c = false;
  }

  if (length == 3) {
    if (all_c) {
      if (s) (*s)->AddNote("{}e - all C penalty (length 3)", hairpin_c3_loop);
      energy += hairpin_c3_loop;
    }
    return energy;
  }
  const Base left = r[st + 1];
  const Base right = r[en - 1];
  if (s) (*s)->AddNote("{}e - terminal mismatch", terminal[r[st]][left][right][r[en]]);
  energy += terminal[r[st]][left][right][r[en]];
  if (IsPairOf(left, right, U_b, U_b) || IsPairOf(left, right, G_b, A_b)) {
    if (s) (*s)->AddNote("{}e - UU/GA first mismatch", hairpin_uu_ga_first_mismatch);
    energy += hairpin_uu_ga_first_mismatch;
  }
  if (IsPairOf(left, right, G_b, G_b)) {
    if (s) (*s)->AddNote("{}e - GG first mismatch", hairpin_gg_first_mismatch);
    energy += hairpin_gg_first_mismatch;
  }
  if (all_c) {
    Energy all_c_energy = hairpin_all_c_a * length + hairpin_all_c_b;
    if (s) (*s)->AddNote("{}e - all C penalty", all_c_energy);
    energy += all_c_energy;
  }
  if (IsPairOf(r[st], r[en], G_b, U_b) && st >= 2 && r[st - 1] == G && r[st - 2] == G) {
    if (s) (*s)->AddNote("{}e - special GU closure", hairpin_special_gu_closure);
    energy += hairpin_special_gu_closure;
  }

  return energy;
}

// Indices are inclusive.
// Rules for bulge loop energy:
// 1. If length (number of unpaired bases in the bulge loop) is more than 1, just use initiation
// energy.
// 2. Otherwise:
//    Don't apply AU/GU penalties -- since the helix is able to continue (unlike for lengths > 1).
//    Since the helix continues, also apply stacking energies for Watson-Crick helices.
//    If the unpaired base is a C, and is next to another C (at pos - 1 or pos + 1), add special C
//    bulge bonus.
//    Count up the number of contiguous bases next to the size 1 bulge loop base, and compute a
//    bonus from that.
Energy Model::Bulge(
    const Primary& r, int ost, int oen, int ist, int ien, std::unique_ptr<Structure>* s) const {
  assert(ist > ost && ien < oen && (oen - ien == 1 || ist - ost == 1) &&
      (oen - ien >= 2 || ist - ost >= 2));
  const int length = std::max(ist - ost, oen - ien) - 1;
  Energy energy = BulgeInitiation(length);

  if (s) {
    *s = std::make_unique<TwoLoopStructure>(ost, oen, ist, ien);
    (*s)->AddNote("Bulge loop len {}", length);
    (*s)->AddNote("{}e - initiation", energy);
  }

  if (length > 1) {
    // Bulges of length > 1 are considered separate helices and get AU/GU penalties.
    if (IsAuPair(r[ost], r[oen])) {
      if (s) (*s)->AddNote("{}e - outer AU penalty", au_penalty);
      energy += au_penalty;
    }
    if (IsGuPair(r[ost], r[oen])) {
      if (s) (*s)->AddNote("{}e - outer GU penalty", gu_penalty);
      energy += gu_penalty;
    }
    if (IsAuPair(r[ist], r[ien])) {
      if (s) (*s)->AddNote("{}e - inner AU penalty", au_penalty);
      energy += au_penalty;
    }
    if (IsGuPair(r[ist], r[ien])) {
      if (s) (*s)->AddNote("{}e - inner GU penalty", gu_penalty);
      energy += gu_penalty;
    }
    return energy;
  }
  // Stacking energy.
  energy += stack[r[ost]][r[ist]][r[ien]][r[oen]];
  if (s) (*s)->AddNote("{}e - stacking", stack[r[ost]][r[ist]][r[ien]][r[oen]]);
  int unpaired = ost + 1;
  if (ost + 1 == ist) unpaired = ien + 1;
  // Special C bulge.
  if (r[unpaired] == C && (r[unpaired - 1] == C || r[unpaired + 1] == C)) {
    if (s) (*s)->AddNote("{}e - special c bulge", bulge_special_c);
    energy += bulge_special_c;
  }

  // Count up the number of contiguous same bases next to the size 1 bulge loop base.
  if (cfg().bulge_states) {
    int num_states = 0;
    for (int i = unpaired; i < static_cast<int>(r.size()) && r[i] == r[unpaired]; ++i) num_states++;
    for (int i = unpaired - 1; i >= 0 && r[i] == r[unpaired]; --i) num_states++;
    Energy states_bonus = -E(R * T * log(num_states));
    if (s) (*s)->AddNote("{}e - {} states bonus", states_bonus, num_states);
    energy += states_bonus;
  }

  return energy;
}

// Indices are inclusive.
// Rules for internal loops:
// 1. If it is 1x1, 1x2, 2x1, or 2x2, then check a lookup table.
// 2. If not, return G_init plus:
// 2.1 Internal loop specific AU/GU penalty for each AU/GU end.
// 2.2 A constant times the absolute difference between the number of unpaired bases on each side of
// the loop.
// 2.3 If the loop is 2x3 or 3x2, look up special mismatch parameters. We just store the values for
// 2x3, and then
//   rotate the rna by 180 degrees to look it up for 3x2.
Energy Model::InternalLoop(
    const Primary& r, int ost, int oen, int ist, int ien, std::unique_ptr<Structure>* s) const {
  const int toplen = ist - ost - 1;
  const int botlen = oen - ien - 1;
  if (s) {
    *s = std::make_unique<TwoLoopStructure>(ost, oen, ist, ien);
    (*s)->AddNote("{}x{} internal loop", toplen, botlen);
  }

  Energy energy = ZERO_E;

  // AU/GU penalties.
  if (IsAuPair(r[ost], r[oen])) {
    if (s) (*s)->AddNote("{}e - outer AU penalty", au_penalty);
    energy += au_penalty;
  }
  if (IsGuPair(r[ost], r[oen])) {
    if (s) (*s)->AddNote("{}e - outer GU penalty", gu_penalty);
    energy += gu_penalty;
  }
  if (IsAuPair(r[ist], r[ien])) {
    if (s) (*s)->AddNote("{}e - inner AU penalty", au_penalty);
    energy += au_penalty;
  }
  if (IsGuPair(r[ist], r[ien])) {
    if (s) (*s)->AddNote("{}e - inner GU penalty", gu_penalty);
    energy += gu_penalty;
  }

  if (toplen == 1 && botlen == 1)
    return energy + internal_1x1[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[oen]];
  if (toplen == 1 && botlen == 2)
    return energy +
        internal_1x2[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]][r[oen]];
  if (toplen == 2 && botlen == 1)
    return energy +
        internal_1x2[r[ien]][r[ien + 1]][r[oen]][r[ost]][r[ost + 1]][r[ost + 2]][r[ist]];
  if (toplen == 2 && botlen == 2)
    return energy +
        internal_2x2[r[ost]][r[ost + 1]][r[ost + 2]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]]
                    [r[oen]];

  // Internal loop extra AU/GU penalties.
  if (IsAuPair(r[ost], r[oen])) {
    if (s) (*s)->AddNote("{}e - outer internal loop AU penalty", internal_au_penalty);
    energy += internal_au_penalty;
  }
  if (IsGuPair(r[ost], r[oen])) {
    if (s) (*s)->AddNote("{}e - outer internal loop GU penalty", internal_gu_penalty);
    energy += internal_gu_penalty;
  }
  if (IsAuPair(r[ist], r[ien])) {
    if (s) (*s)->AddNote("{}e - inner internal loop AU penalty", internal_au_penalty);
    energy += internal_au_penalty;
  }
  if (IsGuPair(r[ist], r[ien])) {
    if (s) (*s)->AddNote("{}e - inner internal loop GU penalty", internal_gu_penalty);
    energy += internal_gu_penalty;
  }

  const auto initiation = InternalLoopInitiation(toplen + botlen);
  energy += initiation;
  if (s) (*s)->AddNote("{}e - initiation", initiation);

  // Asymmetry term, limit with Ninio maximum asymmetry.
  const Energy asym = std::min(std::abs(toplen - botlen) * internal_asym, NINIO_MAX_ASYM);
  if (s) (*s)->AddNote("{}e - asymmetry", asym);
  energy += asym;

  // Special mismatch parameters. To flip an RNA, we flip it vertically and
  // horizontally (180 degree rotation). It turns out that the accesses for 2x3
  // and 3x2 both touch the same array locations, just which side is on the left
  // / right is flipped.
  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2)) {
    const Energy mismatch = internal_2x3_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
        internal_2x3_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];
    if (s) (*s)->AddNote("{}e - 2x3 mismatch params", mismatch);
    energy += mismatch;
  } else if (toplen != 1 && botlen != 1) {
    const Energy mismatch = internal_other_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
        internal_other_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];
    if (s) (*s)->AddNote("{}e - other mismatch params", mismatch);
    energy += mismatch;
  }

  return energy;
}

Energy Model::TwoLoop(
    const Primary& r, int ost, int oen, int ist, int ien, std::unique_ptr<Structure>* s) const {
  const int toplen = ist - ost - 1;
  const int botlen = oen - ien - 1;
  if (toplen == 0 && botlen == 0) {
    if (s) *s = std::make_unique<StackingStructure>(ost, oen);
    return stack[r[ost]][r[ist]][r[ien]][r[oen]];
  }
  if (toplen >= 1 && botlen >= 1) return InternalLoop(r, ost, oen, ist, ien, s);
  return Bulge(r, ost, oen, ist, ien, s);
}

Energy Model::MultiloopEnergy(const Primary& r, const Secondary& s, int st, int en,
    std::deque<int>* branches, bool use_given_ctds, Ctds* ctd,
    std::unique_ptr<Structure>* sstruc) const {
  const bool exterior_loop = s[st] != en;
  Energy energy = ZERO_E;

  std::unique_ptr<MultiLoopStructure> struc = nullptr;
  if (sstruc) {
    struc = std::make_unique<MultiLoopStructure>(st, en);
    if (exterior_loop) struc->AddNote("exterior loop");
  }

  // Add AUGU penalties.
  int num_unpaired = 0;
  for (auto branch_st : *branches) {
    num_unpaired += s[branch_st] - branch_st + 1;

    if (IsAuPair(r[branch_st], r[s[branch_st]])) {
      if (struc)
        struc->AddNote("{}e - opening AU penalty at {} {}", au_penalty, branch_st, s[branch_st]);
      energy += au_penalty;
    }
    if (IsGuPair(r[branch_st], r[s[branch_st]])) {
      if (struc)
        struc->AddNote("{}e - opening GU penalty at {} {}", gu_penalty, branch_st, s[branch_st]);
      energy += gu_penalty;
    }
  }
  num_unpaired = en - st - 1 - num_unpaired + static_cast<int>(exterior_loop) * 2;
  if (struc) struc->AddNote("Unpaired: {}, Branches: {}", num_unpaired, branches->size() + 1);

  BranchCtd branch_ctd;
  Energy ctd_energy = ZERO_E;
  if (exterior_loop) {
    // No initiation for the exterior loop.
    if (use_given_ctds) {
      ctd_energy = AddBaseCtdsToBranchCtds(*this, r, s, *ctd, *branches, &branch_ctd);
    } else {
      ctd_energy = ComputeOptimalCtds(*this, r, s, *branches, true, &branch_ctd);
      AddBranchCtdsToBaseCtds(*branches, branch_ctd, ctd);
    }
  } else {
    if (IsAuPair(r[st], r[en])) {
      if (struc) struc->AddNote("{}e - closing AU penalty at {} {}", au_penalty, st, en);
      energy += au_penalty;
    }
    if (IsGuPair(r[st], r[en])) {
      if (struc) struc->AddNote("{}e - closing GU penalty at {} {}", gu_penalty, st, en);
      energy += gu_penalty;
    }
    Energy initiation = MultiloopInitiation(static_cast<int>(branches->size() + 1));
    if (struc) struc->AddNote("{}e - initiation", initiation);
    energy += initiation;

    if (use_given_ctds) {
      branches->push_front(en);
      ctd_energy = AddBaseCtdsToBranchCtds(*this, r, s, *ctd, *branches, &branch_ctd);
      branches->pop_front();
    } else {
      BranchCtd config_ctds[4] = {};
      std::pair<Energy, int> config_energies[4] = {};
      branches->push_front(en);
      config_energies[0] = {ComputeOptimalCtds(*this, r, s, *branches, true, &config_ctds[0]), 0};
      config_energies[1] = {ComputeOptimalCtds(*this, r, s, *branches, false, &config_ctds[1]), 1};
      branches->pop_front();
      branches->push_back(en);
      config_energies[2] = {ComputeOptimalCtds(*this, r, s, *branches, true, &config_ctds[2]), 2};
      // Swap the final branch back to the front because following code expects it.
      config_ctds[2].push_front(config_ctds[2].back());
      config_ctds[2].pop_back();
      config_energies[3] = {ComputeOptimalCtds(*this, r, s, *branches, false, &config_ctds[3]), 3};
      config_ctds[3].push_front(config_ctds[3].back());
      config_ctds[3].pop_back();
      branches->pop_back();
      std::sort(config_energies, config_energies + 4);
      branch_ctd = config_ctds[config_energies[0].second];
      ctd_energy = config_energies[0].first;

      // Write the optimal ctds to `ctd`.
      branches->push_front(en);
      AddBranchCtdsToBaseCtds(*branches, branch_ctd, ctd);
      branches->pop_front();
    }
  }
  energy += ctd_energy;

  if (struc) {
    struc->AddNote("{}e - ctd", ctd_energy);
    if (!exterior_loop) {
      struc->AddNote("{}e - outer loop - {}", branch_ctd[0].second, CtdToName(branch_ctd[0].first));
      branch_ctd.pop_front();
    }
    for (const auto& [c, e] : branch_ctd) struc->AddCtd(c, e);
    // Give the pointer back.
    *sstruc = std::move(struc);
  }

  return energy;
}

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
// If `ctd` is non-null, use the given ctds.
EnergyResult Model::SubEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd, int st,
    int en, bool build_structure) const {
  ModelBase::VerifyForEfn(r, s, given_ctd);
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

}  // namespace mrna::md::base::opt
