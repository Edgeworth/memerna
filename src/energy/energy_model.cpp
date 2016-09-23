#include "energy/energy_model.h"
#include "energy/structure.h"
#include "parsing.h"

namespace memerna {
namespace energy {

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
energy_t EnergyModel::Hairpin(
    const primary_t& r, int st, int en, std::unique_ptr<Structure>* s) const {
  assert(st < en);
  if (s) *s = std::make_unique<HairpinLoopStructure>(st, en);

  std::string seq;
  for (int i = st; i <= en; ++i) seq += BaseToChar(r[i]);
  const auto iter = hairpin.find(seq);
  if (iter != hairpin.end()) {
    if (s) (*s)->AddNote("special hairpin");
    return iter->second;
  }

  // Subtract two for the initiating base pair.
  const int length = en - st - 1;
  if (length < 3) return MAX_E;  // Disallowed by T04.
  energy_t energy = HairpinInitiation(length);
  if (s) (*s)->AddNote("%de - initiation", energy);
  // Apply AU penalty if necessary (N.B. not for special hairpin sequences).
  if (IsAuGu(r[st], r[en])) {
    if (s) (*s)->AddNote("%de - AU/GU penalty", augu_penalty);
    energy += augu_penalty;
  }

  // T04 says hairpin loops with all C bases inside them are treated specially.
  bool all_c = true;
  for (int i = st + 1; i <= en - 1; ++i) {
    if (r[i] != C) all_c = false;
  }

  if (length == 3) {
    if (all_c) {
      if (s) (*s)->AddNote("%de - all C penalty (length 3)", hairpin_c3_loop);
      energy += hairpin_c3_loop;
    }
    return energy;
  }
  const base_t left = r[st + 1], right = r[en - 1];
  if (s) (*s)->AddNote("%de - terminal mismatch", terminal[r[st]][left][right][r[en]]);
  energy += terminal[r[st]][left][right][r[en]];
  if (IsPairOf(left, right, U_b, U_b) || IsPairOf(left, right, G_b, A_b)) {
    if (s) (*s)->AddNote("%de - UU/GA first mismatch", hairpin_uu_ga_first_mismatch);
    energy += hairpin_uu_ga_first_mismatch;
  }
  if (IsPairOf(left, right, G_b, G_b)) {
    if (s) (*s)->AddNote("%de - GG first mismatch", hairpin_gg_first_mismatch);
    energy += hairpin_gg_first_mismatch;
  }
  if (all_c) {
    energy_t all_c_energy = hairpin_all_c_a * length + hairpin_all_c_b;
    if (s) (*s)->AddNote("%de - all C penalty", all_c_energy);
    energy += all_c_energy;
  }

  if (IsPairOf(r[st], r[en], G_b, U_b) && st >= 2 && r[st - 1] == G && r[st - 2] == G) {
    if (s) (*s)->AddNote("%de - special GU closure", hairpin_special_gu_closure);
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
energy_t EnergyModel::Bulge(
    const primary_t& r, int ost, int oen, int ist, int ien, std::unique_ptr<Structure>* s) const {
  assert(ist > ost && ien < oen && (oen - ien == 1 || ist - ost == 1) &&
         (oen - ien >= 2 || ist - ost >= 2));
  const int length = std::max(ist - ost, oen - ien) - 1;
  energy_t energy = BulgeInitiation(length);

  if (s) {
    *s = std::make_unique<InternalLoopStructure>(ost, oen, ist, ien);
    (*s)->AddNote("Bulge loop len %d", length);
    (*s)->AddNote("%de - initiation", energy);
  }

  if (length > 1) {
    // Bulges of length > 1 are considered separate helices and get AU/GU penalties.
    if (IsAuGu(r[ost], r[oen])) {
      if (s) (*s)->AddNote("%de - outer AU/GU penalty", augu_penalty);
      energy += augu_penalty;
    }
    if (IsAuGu(r[ist], r[ien])) {
      if (s) (*s)->AddNote("%de - inner AU/GU penalty", augu_penalty);
      energy += augu_penalty;
    }
    return energy;
  }
  // Stacking energy.
  energy += stack[r[ost]][r[ist]][r[ien]][r[oen]];
  int unpaired = ost + 1;
  if (ost + 1 == ist) unpaired = ien + 1;
  // Special C bulge.
  if (r[unpaired] == C && (r[unpaired - 1] == C || r[unpaired + 1] == C)) {
    if (s) (*s)->AddNote("%de - special c bulge", bulge_special_c);
    energy += bulge_special_c;
  }

  // Count up the number of contiguous same bases next to the size 1 bulge loop base.
  int num_states = 0;
  for (int i = unpaired; i < int(r.size()) && r[i] == r[unpaired]; ++i) num_states++;
  for (int i = unpaired - 1; i >= 0 && r[i] == r[unpaired]; --i) num_states++;
  energy_t states_bonus = -energy_t(round(10.0 * R * T * log(num_states)));
  if (s) (*s)->AddNote("%de - %d states bonus", states_bonus, num_states);
  energy += states_bonus;

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
energy_t EnergyModel::InternalLoop(
    const primary_t& r, int ost, int oen, int ist, int ien, std::unique_ptr<Structure>* s) const {
  const int toplen = ist - ost - 1, botlen = oen - ien - 1;
  if (s) {
    *s = std::make_unique<InternalLoopStructure>(ost, oen, ist, ien);
    (*s)->AddNote("%dx%d internal loop", toplen, botlen);
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

  energy_t energy = InternalLoopInitiation(toplen + botlen);
  if (s) (*s)->AddNote("%de - initiation", energy);

  // Special AU/GU penalties.
  if (IsAuGu(r[ost], r[oen])) {
    if (s) (*s)->AddNote("%de - outer AU/GU penalty", internal_augu_penalty);
    energy += internal_augu_penalty;
  }
  if (IsAuGu(r[ist], r[ien])) {
    if (s) (*s)->AddNote("%de - inner AU/GU penalty", internal_augu_penalty);
    energy += internal_augu_penalty;
  }
  // Asymmetry term, limit with Ninio maximum asymmetry.
  const energy_t asym = std::min(std::abs(toplen - botlen) * internal_asym, NINIO_MAX_ASYM);
  if (s) (*s)->AddNote("%de - asymmetry", asym);
  energy += asym;

  // Special mismatch parameters.
  // To flip an RNA, we flip it vertically and horizontally (180 degree rotation).
  // It turns out that the accesses for 2x3 and 3x2 both touch the same array locations, just which
  // side is
  // on the left / right is flipped.
  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2)) {
    const energy_t mismatch = internal_2x3_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
                              internal_2x3_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];
    if (s) (*s)->AddNote("%de - 2x3 mismatch params", mismatch);
    energy += mismatch;
  } else if (toplen != 1 && botlen != 1) {
    const energy_t mismatch = internal_other_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
                              internal_other_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];
    if (s) (*s)->AddNote("%de - other mismatch params", mismatch);
    energy += mismatch;
  }

  return energy;
}

energy_t EnergyModel::TwoLoop(
    const primary_t& r, int ost, int oen, int ist, int ien, std::unique_ptr<Structure>* s) const {
  const int toplen = ist - ost - 1, botlen = oen - ien - 1;
  if (toplen == 0 && botlen == 0) {
    if (s) *s = std::make_unique<StackingStructure>(ost, oen);
    return stack[r[ost]][r[ist]][r[ien]][r[oen]];
  }
  if (toplen >= 1 && botlen >= 1) return InternalLoop(r, ost, oen, ist, ien, s);
  return Bulge(r, ost, oen, ist, ien, s);
}

uint32_t EnergyModel::Checksum() const {
  std::string data;

// This isn't portable across machines with different endianness but I don't care.
#define APPEND_DATA(d)                           \
  do {                                           \
    auto dp = reinterpret_cast<const char*>(&d); \
    data.insert(data.end(), dp, dp + sizeof(d)); \
  } while (0)

  APPEND_DATA(stack);
  APPEND_DATA(terminal);
  APPEND_DATA(internal_init);
  APPEND_DATA(internal_1x1);
  APPEND_DATA(internal_1x2);
  APPEND_DATA(internal_2x2);
  APPEND_DATA(internal_2x3_mismatch);
  APPEND_DATA(internal_other_mismatch);
  APPEND_DATA(internal_asym);
  APPEND_DATA(internal_augu_penalty);
  APPEND_DATA(bulge_init);
  APPEND_DATA(bulge_special_c);
  APPEND_DATA(hairpin_init);
  APPEND_DATA(hairpin_uu_ga_first_mismatch);
  APPEND_DATA(hairpin_gg_first_mismatch);
  APPEND_DATA(hairpin_special_gu_closure);
  APPEND_DATA(hairpin_c3_loop);
  APPEND_DATA(hairpin_all_c_a);
  APPEND_DATA(hairpin_all_c_b);

  for (const auto& v : hairpin) {
    data += v.first;
    APPEND_DATA(v.second);
  }

  APPEND_DATA(multiloop_hack_a);
  APPEND_DATA(multiloop_hack_b);
  APPEND_DATA(dangle5);
  APPEND_DATA(dangle3);
  APPEND_DATA(coax_mismatch_non_contiguous);
  APPEND_DATA(coax_mismatch_wc_bonus);
  APPEND_DATA(coax_mismatch_gu_bonus);
  APPEND_DATA(augu_penalty);

  APPEND_DATA(HAIRPIN_MIN_SZ);
  APPEND_DATA(R);
  APPEND_DATA(T);
  APPEND_DATA(NINIO_MAX_ASYM);
  APPEND_DATA(TWOLOOP_MAX_SZ);
#undef APPEND_DATA

  return Crc32(data);
}

bool EnergyModel::IsValid(std::string* reason) const {
#define CHECK_COND(cond, reason_str)                           \
  do {                                                         \
    if (!(cond)) {                                             \
      if (reason) *reason = "expected " #cond ": " reason_str; \
      return false;                                            \
    }                                                          \
  } while (0)

  for (int a = 0; a < 4; ++a) {
    for (int b = 0; b < 4; ++b) {
      for (int c = 0; c < 4; ++c) {
        for (int d = 0; d < 4; ++d) {
          // Expect 180 degree rotations to be the same.
          CHECK_COND(
              stack[a][b][c][d] == stack[c][d][a][b], "180 degree rotations should be the same");
          CHECK_COND(internal_asym >= 0, "optimisations rely on this");

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
#undef CHECK_COND
  return true;
}
}
}
