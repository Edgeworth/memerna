#include <cstdio>
#include <cmath>
#include <memory>
#include <algorithm>
#include "energy/energy.h"
#include "energy/structure.h"

namespace memerna {
namespace energy {

using namespace internal;

namespace {
std::vector<Ctd> g_ctds;
}

// Indices are inclusive, include the initiating base pair.
// N.B. This includes an ending AU/GU penalty.
// Rules for hairpin energy:
// 1. Check a lookup table for a hardcoded value.
// 2. Loops with less than 3 bases on the inside are disallowed by T04.
// Case 1 -- 3 bases on the inside:
//   Return G_init + penalty if all bases inside are C.
//   The penalty for all C bases for size 3 hairpin loops is special cased.
// Case 2 -- more than 3 bases on the inside:
//   Return G_init plus the following penalties / bonuses:
//   Terminal mismatch energy for st + 1 and en - 1.
//   If the mismatch is UU or GA (not AG), additional bonus
//   If the mismatch is GG, additional bonus.
//   If the pair st, en is GU (not UG), a bonus if st - 1 and st - 2 are both Gs, if they exist.
//   A penalty if all the bases inside are C: A * length + B (A, B specified as part of the energy model).
energy_t Hairpin(const primary_t& r, int st, int en, std::unique_ptr<Structure>* s) {
  assert(st < en);
  if (s) *s = std::make_unique<HairpinLoopStructure>(st, en);

  std::string seq;
  for (int i = st; i <= en; ++i)
    seq += BaseToChar(r[i]);
  auto iter = g_hairpin.find(seq);
  if (iter != g_hairpin.end()) {
    if (s) (*s)->AddNote("special hairpin");
    return iter->second;
  }

  // Subtract two for the initiating base pair.
  int length = en - st - 1;
  if (length < 3) return constants::MAX_E;  // Disallowed by T04.
  energy_t energy = HairpinInitiation(length);
  if (s) (*s)->AddNote("%de - initiation", energy);
  // Apply AU penalty if necessary (N.B. not for special hairpin sequences).
  if (IsAuGu(r[st], r[en])) {
    if (s) (*s)->AddNote("%de - AU/GU penalty", g_augu_penalty);
    energy += g_augu_penalty;
  }

  // T04 says hairpin loops with all C bases inside them are treated specially.
  bool all_c = true;
  for (int i = st + 1; i <= en - 1; ++i) {
    if (r[i] != C) all_c = false;
  }

  if (length == 3) {
    if (all_c) {
      if (s) (*s)->AddNote("%de - all C penalty (length 3)", g_hairpin_c3_loop);
      energy += g_hairpin_c3_loop;
    }
    return energy;
  }
  base_t left = r[st + 1], right = r[en - 1];
  if (s) (*s)->AddNote("%de - terminal mismatch", g_terminal[r[st]][left][right][r[en]]);
  energy += g_terminal[r[st]][left][right][r[en]];
  if (IsPairOf(left, right, U_b, U_b) || IsPairOf(left, right, G_b, A_b)) {
    if (s) (*s)->AddNote("%de - UU/GA first mismatch", g_hairpin_uu_ga_first_mismatch);
    energy += g_hairpin_uu_ga_first_mismatch;
  }
  if (IsPairOf(left, right, G_b, G_b)) {
    if (s) (*s)->AddNote("%de - GG first mismatch", g_hairpin_gg_first_mismatch);
    energy += g_hairpin_gg_first_mismatch;
  }
  if (all_c) {
    energy_t all_c_energy = g_hairpin_all_c_a * length + g_hairpin_all_c_b;
    if (s) (*s)->AddNote("%de - all C penalty", all_c_energy);
    energy += all_c_energy;
  }

  if (IsPairOf(r[st], r[en], G_b, U_b) && st >= 2 && r[st - 1] == G && r[st - 2] == G) {
    if (s) (*s)->AddNote("%de - special GU closure", g_hairpin_special_gu_closure);
    energy += g_hairpin_special_gu_closure;
  }
  return energy;
}

// Indices are inclusive.
// Rules for bulge loop energy:
// 1. If length (number of unpaired bases in the bulge loop) is more than 1, just use initiation energy.
// 2. Otherwise:
//    Don't apply AU/GU penalties -- since the helix is able to continue (unlike for lengths > 1).
//    Since the helix continues, also apply stacking energies for Watson-Crick helices.
//    If the unpaired base is a C, and is next to another C (at pos - 1 or pos + 1), add special C bulge bonus.
//    Count up the number of contiguous bases next to the size 1 bulge loop base, and compute a bonus from that.
energy_t Bulge(const primary_t& r,
    int ost, int oen, int ist, int ien, std::unique_ptr<Structure>* s) {
  assert(ist > ost && ien < oen && (oen - ien == 1 || ist - ost == 1) && (oen - ien >= 2 || ist - ost >= 2));
  int length = std::max(ist - ost, oen - ien) - 1;
  energy_t energy = BulgeInitiation(length);

  if (s) {
    *s = std::make_unique<InternalLoopStructure>(ost, oen, ist, ien);
    (*s)->AddNote("Bulge loop len %d", length);
    (*s)->AddNote("%de - initiation", energy);
  }

  if (length > 1) {
    // Bulges of length > 1 are considered separate helices and get AU/GU penalties.
    if (IsAuGu(r[ost], r[oen])) {
      if (s) (*s)->AddNote("%de - outer AU/GU penalty", g_augu_penalty);
      energy += g_augu_penalty;
    }
    if (IsAuGu(r[ist], r[ien])) {
      if (s) (*s)->AddNote("%de - inner AU/GU penalty", g_augu_penalty);
      energy += g_augu_penalty;
    }
    return energy;
  }
  // Stacking energy.
  energy += g_stack[r[ost]][r[ist]][r[ien]][r[oen]];
  int unpaired = ost + 1;
  if (ost + 1 == ist) unpaired = ien + 1;
  // Special C bulge.
  if (r[unpaired] == C && (r[unpaired - 1] == C || r[unpaired + 1] == C)) {
    if (s) (*s)->AddNote("%de - special c bulge", g_bulge_special_c);
    energy += g_bulge_special_c;
  }

  // Count up the number of contiguous same bases next to the size 1 bulge loop base.
  int num_states = 0;
  for (int i = unpaired; i < int(r.size()) && r[i] == r[unpaired]; ++i)
    num_states++;
  for (int i = unpaired - 1; i >= 0 && r[i] == r[unpaired]; --i)
    num_states++;
  energy_t states_bonus = -energy_t(round(10.0 * constants::R * constants::T * log(num_states)));
  if (s) (*s)->AddNote("%de - %d states bonus", states_bonus, num_states);
  energy += states_bonus;

  return energy;
}

// Indices are inclusive.
// Rules for internal loops:
// 1. If it is 1x1, 1x2, 2x1, or 2x2, then check a lookup table.
// 2. If not, return G_init plus:
// 2.1 Internal loop specific AU/GU penalty for each AU/GU end.
// 2.2 A constant times the absolute difference between the number of unpaired bases on each side of the loop.
// 2.3 If the loop is 2x3 or 3x2, look up special mismatch parameters. We just store the values for 2x3, and then
//   rotate the rna by 180 degrees to look it up for 3x2.
energy_t InternalLoop(const primary_t& r,
    int ost, int oen, int ist, int ien, std::unique_ptr<Structure>* s) {
  int toplen = ist - ost - 1, botlen = oen - ien - 1;
  if (s) {
    *s = std::make_unique<InternalLoopStructure>(ost, oen, ist, ien);
    (*s)->AddNote("%dx%d internal loop", toplen, botlen);
  }
  if (toplen == 1 && botlen == 1)
    return g_internal_1x1[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[oen]];
  if (toplen == 1 && botlen == 2)
    return g_internal_1x2[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]][r[oen]];
  if (toplen == 2 && botlen == 1)
    return g_internal_1x2[r[ien]][r[ien + 1]][r[oen]][r[ost]][r[ost + 1]][r[ost + 2]][r[ist]];
  if (toplen == 2 && botlen == 2)
    return g_internal_2x2[r[ost]][r[ost + 1]][r[ost + 2]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]][r[oen]];

  energy_t energy = InternalLoopInitiation(toplen + botlen);
  if (s) (*s)->AddNote("%de - initiation", energy);

  // Special AU/GU penalties.
  if (IsAuGu(r[ost], r[oen])) {
    if (s) (*s)->AddNote("%de - outer AU/GU penalty", g_internal_augu_penalty);
    energy += g_internal_augu_penalty;
  }
  if (IsAuGu(r[ist], r[ien])) {
    if (s) (*s)->AddNote("%de - inner AU/GU penalty", g_internal_augu_penalty);
    energy += g_internal_augu_penalty;
  }
  // Asymmetry term, limit with Ninio maximum asymmetry.
  energy_t asym = std::min(std::abs(toplen - botlen) * g_internal_asym, constants::NINIO_MAX_ASYM);
  if (s) (*s)->AddNote("%de - asymmetry", asym);
  energy += asym;

  // Special mismatch parameters.
  // To flip an RNA, we flip it vertically and horizontally (180 degree rotation).
  // It turns out that the accesses for 2x3 and 3x2 both touch the same array locations, just which side is
  // on the left / right is flipped.
  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2)) {
    energy_t mismatch =
        g_internal_2x3_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
            g_internal_2x3_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];
    if (s) (*s)->AddNote("%de - 2x3 mismatch params", mismatch);
    energy += mismatch;
  } else if (toplen != 1 && botlen != 1) {
    energy_t mismatch =
        g_internal_other_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
            g_internal_other_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];
    if (s) (*s)->AddNote("%de - other mismatch params", mismatch);
    energy += mismatch;
  }

  return energy;
}

energy_t TwoLoop(const primary_t& r,
    int ost, int oen, int ist, int ien, std::unique_ptr<Structure>* s) {
  int toplen = ist - ost - 1, botlen = oen - ien - 1;
  if (toplen == 0 && botlen == 0) {
    if (s) *s = std::make_unique<StackingStructure>(ost, oen);
    return g_stack[r[ost]][r[ist]][r[ien]][r[oen]];
  }
  if (toplen >= 1 && botlen >= 1)
    return InternalLoop(r, ost, oen, ist, ien, s);
  return Bulge(r, ost, oen, ist, ien, s);
}

energy_t MultiloopEnergy(const secondary_t& secondary,
    int st, int en, std::deque<int>& branches, std::unique_ptr<Structure>* s) {
  const auto& r = secondary.r;
  const auto& p = secondary.p;
  bool exterior_loop = st == 0 && en == int(r.size() - 1) && p[st] != en;
  energy_t energy = 0;

  if (s) {
    *s = std::make_unique<MultiLoopStructure>(st, en);
    if (exterior_loop) (*s)->AddNote("exterior loop");
  }

  // Add AUGU penalties.
  int num_unpaired = 0;
  for (auto branch_st : branches) {
    num_unpaired += p[branch_st] - branch_st + 1;

    if (IsAuGu(r[branch_st], r[p[branch_st]])) {
      if (s) (*s)->AddNote("%de - opening AU/GU penalty at %d %d", g_augu_penalty, branch_st, p[branch_st]);
      energy += g_augu_penalty;
    }
  }
  num_unpaired = en - st - 1 - num_unpaired + exterior_loop * 2;
  if (s) (*s)->AddNote("Unpaired: %d, Branches: %zu", num_unpaired, branches.size() + 1);

  bool compute_ctd = branches.empty() || (g_ctds[branches.front()] == CTD_NA);
  branch_ctd_t branch_ctds;
  energy_t ctd_energy = 0;
  if (exterior_loop) {
    // No initiation for the exterior loop.
    if (compute_ctd) {
      ctd_energy = ComputeOptimalCtd(secondary, branches, -1, true, &branch_ctds);
      BranchCtdsToBaseCtds(secondary, branches, branch_ctds, &g_ctds);
    } else {
      BaseCtdsToBranchCtds(secondary, branches, g_ctds, &branch_ctds);
    }
  } else {
    if (IsAuGu(r[st], r[en])) {
      if (s) (*s)->AddNote("%de - closing AU/GU penalty at %d %d", g_augu_penalty, st, en);
      energy += g_augu_penalty;
    }
    energy_t initiation = MultiloopInitiation(int(branches.size() + 1));
    if (s) (*s)->AddNote("%de - initiation", initiation);
    energy += initiation;

    if (compute_ctd) {
      branch_ctd_t config_ctds[4] = {};
      std::pair<energy_t, int> config_energies[4] = {};
      branches.push_front(st);
      config_energies[0] = {ComputeOptimalCtd(secondary, branches, 0, true, config_ctds), 0};
      config_energies[1] = {ComputeOptimalCtd(secondary, branches, 0, false, config_ctds + 1), 1};
      branches.pop_front();
      branches.push_back(st);
      config_energies[2] = {ComputeOptimalCtd(
          secondary, branches, int(branches.size() - 1), true, config_ctds + 2), 2};
      config_energies[3] = {ComputeOptimalCtd(
          secondary, branches, int(branches.size() - 1), false, config_ctds + 3), 3};
      // Shuffle these around so the outer loop is as the start.
      config_ctds[2].push_front(config_ctds[2].back());
      config_ctds[2].pop_back();
      config_ctds[3].push_front(config_ctds[3].back());
      config_ctds[3].pop_back();
      branches.pop_back();
      std::sort(config_energies, config_energies + 4);
      branch_ctds = config_ctds[config_energies[0].second];
      ctd_energy = config_energies[0].first;

      // Assumes the outer loop is the first one.
      // TODO uncomment.
//      branches.push_front(st);
//      BranchCtdsToBaseCtds(branches, branch_ctds);
//#ifndef NDEBUG
//      branch_ctd_t tmp;
//      assert(ctd_energy == BaseCtdsToBranchCtds(branches, &tmp));
//      assert(branch_ctds == tmp);
//#endif
//      branches.pop_front();
    } else {
      branches.push_front(st);
      BaseCtdsToBranchCtds(secondary, branches, g_ctds, &branch_ctds);
      branches.pop_front();
    }
  }
  energy += ctd_energy;

  if (s) {
    (*s)->AddNote("%de - ctd", ctd_energy);
    // TODO uncomment.
//    int idx = 0;
//    if (!exterior_loop) {
//      (*s)->AddNote("%de - outer loop stacking - %s", branch_ctds[0].second,
//          energy::CtdToName(branch_ctds[0].first));
//      idx++;
//    }
//    for (; idx < branch_ctds.size(); ++i) {
//      (*s)->AddBranch()
//    }
    // TODO note for outer loop stacking
  }

  return energy;
}

energy_t ComputeEnergyInternal(const secondary_t& secondary, int st, int en, std::unique_ptr<Structure>* s) {
  const auto& r = secondary.r;
  const auto& p = secondary.p;
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
  bool exterior_loop = st == 0 && en == int(r.size() - 1) && p[st] != en;
  if (exterior_loop || branches.size() >= 2) {
    // Multiloop.
    energy += MultiloopEnergy(secondary, st, en, branches, s);  // TODO add branch occurs below here.
  } else if (branches.size() == 0) {
    // Hairpin loop.
    assert(en - st - 1 >= 3);
    energy += Hairpin(r, st, en, s);
  } else if (branches.size() == 1) {
    int loop_st = branches.front(), loop_en = p[branches.front()];
    energy += TwoLoop(r, st, en, loop_st, loop_en, s);
  }

  if (s) (*s)->SetSelfEnergy(energy);
  // Add energy from children.
  for (auto i : branches) {
    int pair = p[i];
    if (s) {
      std::unique_ptr<Structure> structure;
      energy += ComputeEnergyInternal(secondary, i, pair, &structure);
      (*s)->AddBranch(std::move(structure));
    } else {
      energy += ComputeEnergyInternal(secondary, i, pair, nullptr);
    }
  }
  if (s) (*s)->SetTotalEnergy(energy);

  return energy;
}

computed_t ComputeEnergy(const secondary_t& secondary, std::unique_ptr<Structure>* s) {
  const auto& r = secondary.r;
  const auto& p = secondary.p;
  g_ctds.clear();
  g_ctds.resize(r.size(), CTD_NA);
  assert(r.size() == p.size() && r.size() == g_ctds.size());
  energy_t energy = ComputeEnergyInternal(secondary, 0, (int) r.size() - 1, s);
  if (p[0] == int(r.size() - 1) && IsAuGu(r[0], r[p[0]])) {
    energy += g_augu_penalty;
    if (s) {
      (*s)->AddNote("%de - top level AU/GU penalty", g_augu_penalty);
      (*s)->SetSelfEnergy((*s)->GetSelfEnergy() + g_augu_penalty);  // Gross.
      (*s)->SetTotalEnergy((*s)->GetTotalEnergy() + g_augu_penalty);  // Gross.
    }
  }
  return {r, p, g_ctds, energy};
}

}
}
