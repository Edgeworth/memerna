// Copyright 2023 Eliot Courtney.
#include "compute/energy/common/t04like/model_mixin.h"

#include <fmt/core.h>

#include <algorithm>
#include <compare>
#include <istream>
#include <optional>
#include <utility>

#include "compute/energy/common/branch.h"
#include "compute/energy/common/model.h"
#include "compute/energy/common/parse.h"
#include "compute/energy/common/t04like/branch.h"
#include "compute/energy/energy_cfg.h"
#include "compute/energy/structure.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "util/error.h"
#include "util/string.h"

namespace mrna::erg {

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
//   A penalty if all the bases inside are C: A * length + B (A, B specified as part of the energy
//   model).
Energy T04ModelMixin::Hairpin(
    const Primary& r, int st, int en, std::unique_ptr<Structure>* s) const {
  assert(st < en);
  if (s) *s = std::make_unique<HairpinLoopStructure>(st, en);

  std::string seq;
  for (int i = st; i <= en; ++i) seq += BaseToChar(r[i]).value();
  const auto iter = hairpin.find(seq);
  if (iter != hairpin.end()) {
    if (s) (*s)->AddNote("special hairpin");
    return iter->second;
  }

  // Subtract two for the initiating base pair.
  const int length = en - st - 1;
  if (length < 3) return MAX_E;  // Disallowed by T04.
  Energy energy = HairpinInitiation(length);
  if (s) (*s)->AddNote("{}e - initiation", energy);
  // Apply AU penalty if necessary (N.B. not for special hairpin sequences).
  if (IsAuPair(r[st], r[en])) {
    if (s) (*s)->AddNote("{}e - AU penalty", au_penalty);
    energy += au_penalty;
  }
  if (IsGuPair(r[st], r[en])) {
    if (s) (*s)->AddNote("{}e - GU penalty", gu_penalty);
    energy += gu_penalty;
  }

  // T04 says hairpin loops with all C bases inside them are treated specially.
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
Energy T04ModelMixin::Bulge(
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
  if (cfg.bulge_states) {
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
Energy T04ModelMixin::InternalLoop(
    const Primary& r, int ost, int oen, int ist, int ien, std::unique_ptr<Structure>* s) const {
  const int toplen = ist - ost - 1;
  const int botlen = oen - ien - 1;
  if (s) {
    *s = std::make_unique<TwoLoopStructure>(ost, oen, ist, ien);
    (*s)->AddNote("{}x{} internal loop", toplen, botlen);
  }

  if (toplen == 1 && botlen == 1)
    return internal_1x1[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[oen]];
  if (toplen == 1 && botlen == 2)
    return internal_1x2[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]][r[oen]];
  if (toplen == 2 && botlen == 1)
    return internal_1x2[r[ien]][r[ien + 1]][r[oen]][r[ost]][r[ost + 1]][r[ost + 2]][r[ist]];
  if (toplen == 2 && botlen == 2)
    return internal_2x2[r[ost]][r[ost + 1]][r[ost + 2]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]]
                       [r[oen]];

  Energy energy = InternalLoopInitiation(toplen + botlen);
  if (s) (*s)->AddNote("{}e - initiation", energy);

  // Special AU/GU penalties.
  if (IsAuPair(r[ost], r[oen])) {
    if (s) (*s)->AddNote("{}e - outer AU penalty", internal_au_penalty);
    energy += internal_au_penalty;
  }
  if (IsGuPair(r[ost], r[oen])) {
    if (s) (*s)->AddNote("{}e - outer GU penalty", internal_gu_penalty);
    energy += internal_gu_penalty;
  }
  if (IsAuPair(r[ist], r[ien])) {
    if (s) (*s)->AddNote("{}e - inner AU penalty", internal_au_penalty);
    energy += internal_au_penalty;
  }
  if (IsGuPair(r[ist], r[ien])) {
    if (s) (*s)->AddNote("{}e - inner GU penalty", internal_gu_penalty);
    energy += internal_gu_penalty;
  }
  // Asymmetry term, limit with Ninio maximum asymmetry.
  const Energy asym = std::min(std::abs(toplen - botlen) * internal_asym, NINIO_MAX_ASYM);
  if (s) (*s)->AddNote("{}e - asymmetry", asym);
  energy += asym;

  // Special mismatch parameters.
  // To flip an RNA, we flip it vertically and horizontally (180 degree rotation).
  // It turns out that the accesses for 2x3 and 3x2 both touch the same array locations, just which
  // side is
  // on the left / right is flipped.
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

Energy T04ModelMixin::TwoLoop(
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

Energy T04ModelMixin::MultiloopEnergy(const Primary& r, const Secondary& s, int st, int en,
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

      // Write the optimal ctds to |ctd|.
      branches->push_front(en);
      AddBranchCtdsToBaseCtds(*branches, branch_ctd, ctd);
      branches->pop_front();
    }
  }
  energy += ctd_energy;

  if (struc) {
    struc->AddNote("{}e - ctd", ctd_energy);
    if (!exterior_loop) {
      struc->AddNote(
          "{}e - outer loop - {}", branch_ctd[0].second, erg::CtdToName(branch_ctd[0].first));
      branch_ctd.pop_front();
    }
    for (const auto& [c, e] : branch_ctd) struc->AddCtd(c, e);
    // Give the pointer back.
    *sstruc = std::move(struc);
  }

  return energy;
}

bool T04ModelMixin::IsValid(std::string* reason) const {
  for (int a = 0; a < 4; ++a) {
    for (int b = 0; b < 4; ++b) {
      for (int c = 0; c < 4; ++c) {
        for (int d = 0; d < 4; ++d) {
          // Expect 180 degree rotations to be the same.
          CHECK_COND(
              stack[a][b][c][d] == stack[c][d][a][b], "180 degree rotations should be the same");
          CHECK_COND(internal_asym >= ZERO_E, "optimisations rely on this");

          for (int e = 0; e < 4; ++e) {
            for (int f = 0; f < 4; ++f) {
              CHECK_COND(internal_1x1[a][b][c][d][e][f] == internal_1x1[d][e][f][a][b][c],
                  "180 degree rotations should be the same");
              for (int g = 0; g < 4; ++g) {
                for (int h = 0; h < 4; ++h) {
                  CHECK_COND(
                      internal_2x2[a][b][c][d][e][f][g][h] == internal_2x2[e][f][g][h][a][b][c][d],
                      "180 degree rotations should be the same");
                }
              }
            }
          }
        }
      }
    }
  }
  return true;
}

void T04ModelMixin::ParseMiscDataFromFile(const std::string& filename) {
  std::ifstream f(filename);
  verify(f, "could not open file");

#define READ_DATA(var)                              \
  do {                                              \
    while (1) {                                     \
      auto line = sgetline(f);                      \
      verify(f, "unexpected EOF or error");         \
      if (line.empty() || line[0] == '/') continue; \
      (var) = Energy::FromString(Trim(line));       \
      break;                                        \
    }                                               \
  } while (0)

  // Bulge loops.
  READ_DATA(bulge_special_c);

  // Coaxial stacking.
  READ_DATA(coax_mismatch_non_contiguous);
  READ_DATA(coax_mismatch_wc_bonus);
  READ_DATA(coax_mismatch_gu_bonus);

  // Hairpin loops.
  READ_DATA(hairpin_uu_ga_first_mismatch);
  READ_DATA(hairpin_gg_first_mismatch);
  READ_DATA(hairpin_special_gu_closure);
  READ_DATA(hairpin_c3_loop);
  READ_DATA(hairpin_all_c_a);
  READ_DATA(hairpin_all_c_b);

  // Internal loops.
  READ_DATA(internal_asym);
  READ_DATA(internal_au_penalty);
  READ_DATA(internal_gu_penalty);

  // Multiloop data.
  READ_DATA(multiloop_hack_a);
  READ_DATA(multiloop_hack_b);

  // AU/GU penalty
  READ_DATA(au_penalty);
  READ_DATA(gu_penalty);
#undef READ_DATA

  f >> std::ws;
  verify(f.eof(), "expected EOF");
}

void T04ModelMixin::LoadFromDir(const std::string& data_dir) {
  // Stacking interaction data.
  Parse4MapFromFile(data_dir + "/stacking.data", stack);

  // Terminal mismatch data.
  Parse4MapFromFile(data_dir + "/terminal.data", terminal);

  // Hairpin data.
  ParseNMapFromFile(data_dir + "/hairpin.data", &hairpin);
  ParseVecFromFile(data_dir + "/hairpin_initiation.data", hairpin_init);

  // Bulge loop data.
  ParseVecFromFile(data_dir + "/bulge_initiation.data", bulge_init);

  // Internal loop data.
  ParseVecFromFile(data_dir + "/internal_initiation.data", internal_init);
  Parse6MapFromFile(data_dir + "/internal_1x1.data", internal_1x1);
  Parse7MapFromFile(data_dir + "/internal_1x2.data", internal_1x2);
  Parse8MapFromFile(data_dir + "/internal_2x2.data", internal_2x2);
  Parse4MapFromFile(data_dir + "/internal_2x3_mismatch.data", internal_2x3_mismatch);
  Parse4MapFromFile(data_dir + "/internal_other_mismatch.data", internal_other_mismatch);

  // Dangle data.
  Parse3MapFromFile(data_dir + "/dangle3.data", dangle3);
  Parse3MapFromFile(data_dir + "/dangle5.data", dangle5);

  // Other misc data.
  ParseMiscDataFromFile(data_dir + "/misc.data");
}

void T04ModelMixin::LoadRandom(std::mt19937& eng) {
  std::uniform_real_distribution<double> energy_dist(RAND_MIN_ENERGY, RAND_MAX_ENERGY);
  std::uniform_real_distribution<double> nonneg_energy_dist(0, RAND_MAX_ENERGY);

  RANDOMISE_DATA(stack);
  RANDOMISE_DATA(terminal);
  RANDOMISE_DATA(internal_init);
  RANDOMISE_DATA(internal_1x1);
  RANDOMISE_DATA(internal_1x2);
  RANDOMISE_DATA(internal_2x2);
  RANDOMISE_DATA(internal_2x3_mismatch);
  RANDOMISE_DATA(internal_other_mismatch);
  // This needs to be non-negative for some optimisations.
  internal_asym = Energy::FromDouble(nonneg_energy_dist(eng));
  RANDOMISE_DATA(internal_au_penalty);
  RANDOMISE_DATA(internal_gu_penalty);
  RANDOMISE_DATA(bulge_init);
  RANDOMISE_DATA(bulge_special_c);
  RANDOMISE_DATA(hairpin_init);
  RANDOMISE_DATA(hairpin_uu_ga_first_mismatch);
  RANDOMISE_DATA(hairpin_gg_first_mismatch);
  RANDOMISE_DATA(hairpin_special_gu_closure);
  RANDOMISE_DATA(hairpin_c3_loop);
  RANDOMISE_DATA(hairpin_all_c_a);
  RANDOMISE_DATA(hairpin_all_c_b);

  hairpin.clear();
  std::uniform_int_distribution<int> hairpin_size_dist(HAIRPIN_MIN_SZ, RAND_MAX_HAIRPIN_SZ);
  static_assert(HAIRPIN_MIN_SZ <= RAND_MAX_HAIRPIN_SZ,
      "HAIRPIN_MIN_SZ > RAND_MAX_HAIRPIN does not make sense");
  std::uniform_int_distribution<int> num_hairpin_dist(0, RAND_MAX_NUM_HAIRPIN);
  const int num_hairpin = num_hairpin_dist(eng);
  for (int i = 0; i < num_hairpin; ++i) {
    auto hp = Primary::Random(hairpin_size_dist(eng)).ToSeq();
    hairpin[hp] = Energy::FromDouble(energy_dist(eng));
  }

  RANDOMISE_DATA(multiloop_hack_a);
  RANDOMISE_DATA(multiloop_hack_b);
  RANDOMISE_DATA(dangle5);
  RANDOMISE_DATA(dangle3);
  RANDOMISE_DATA(coax_mismatch_non_contiguous);
  RANDOMISE_DATA(coax_mismatch_wc_bonus);
  RANDOMISE_DATA(coax_mismatch_gu_bonus);
  RANDOMISE_DATA(au_penalty);
  RANDOMISE_DATA(gu_penalty);

  for (int a = 0; a < 4; ++a) {
    for (int b = 0; b < 4; ++b) {
      for (int c = 0; c < 4; ++c) {
        for (int d = 0; d < 4; ++d) {
          // Correct things to be 180 degree rotations if required.
          stack[c][d][a][b] = stack[a][b][c][d];
          for (int e = 0; e < 4; ++e) {
            for (int f = 0; f < 4; ++f) {
              internal_1x1[d][e][f][a][b][c] = internal_1x1[a][b][c][d][e][f];
              for (int g = 0; g < 4; ++g) {
                for (int h = 0; h < 4; ++h) {
                  internal_2x2[e][f][g][h][a][b][c][d] = internal_2x2[a][b][c][d][e][f][g][h];
                }
              }
            }
          }
        }
      }
    }
  }
}

}  // namespace mrna::erg
