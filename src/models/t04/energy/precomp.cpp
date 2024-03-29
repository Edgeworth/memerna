// Copyright 2016 E.
#include "models/t04/energy/precomp.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <compare>
#include <cstdlib>
#include <utility>

#include "model/constants.h"
#include "model/energy.h"
#include "models/common/branch.h"

namespace mrna::md::t04 {

namespace {

Energy MinEnergy(const Energy* energy, std::size_t size) {
  Energy min = energy[0];
  for (int i = 0; i < static_cast<int>(size / sizeof(Energy)); ++i) min = std::min(min, energy[i]);
  return min;
}

}  // namespace

Precomp::Precomp(Primary r, Model::Ptr em) : r_(std::move(r)), em_(std::move(em)) {
  PrecomputeData();
}

Energy Precomp::TwoLoop(int ost, int oen, int ist, int ien) const {
  const int toplen = ist - ost - 1;
  const int botlen = oen - ien - 1;
  if (toplen == 0 && botlen == 0) return em().stack[r_[ost]][r_[ist]][r_[ien]][r_[oen]];
  if (toplen == 0 || botlen == 0) return em().Bulge(r_, ost, oen, ist, ien);

  Energy energy = ZERO_E;
  energy += em().AuGuPenalty(r_[ost], r_[oen]);
  energy += em().AuGuPenalty(r_[ist], r_[ien]);

  if (toplen == 1 && botlen == 1)
    return energy + em().internal_1x1[r_[ost]][r_[ost + 1]][r_[ist]][r_[ien]][r_[ien + 1]][r_[oen]];
  if (toplen == 1 && botlen == 2)
    return energy +
        em().internal_1x2[r_[ost]][r_[ost + 1]][r_[ist]][r_[ien]][r_[ien + 1]][r_[ien + 2]]
                         [r_[oen]];
  if (toplen == 2 && botlen == 1)
    return energy +
        em().internal_1x2[r_[ien]][r_[ien + 1]][r_[oen]][r_[ost]][r_[ost + 1]][r_[ost + 2]]
                         [r_[ist]];
  if (toplen == 2 && botlen == 2)
    return energy +
        em().internal_2x2[r_[ost]][r_[ost + 1]][r_[ost + 2]][r_[ist]][r_[ien]][r_[ien + 1]]
                         [r_[ien + 2]][r_[oen]];

  static_assert(TWOLOOP_MAX_SZ <= Model::INITIATION_CACHE_SZ, "initiation cache not large enough");
  assert(toplen + botlen < Model::INITIATION_CACHE_SZ);
  energy += em().internal_init[toplen + botlen] +
      std::min(std::abs(toplen - botlen) * em().internal_asym, NINIO_MAX_ASYM);

  energy += em().InternalLoopAuGuPenalty(r_[ost], r_[oen]);
  energy += em().InternalLoopAuGuPenalty(r_[ist], r_[ien]);

  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2))
    energy += em().internal_2x3_mismatch[r_[ost]][r_[ost + 1]][r_[oen - 1]][r_[oen]] +
        em().internal_2x3_mismatch[r_[ien]][r_[ien + 1]][r_[ist - 1]][r_[ist]];
  else if (toplen != 1 && botlen != 1)
    energy += em().internal_other_mismatch[r_[ost]][r_[ost + 1]][r_[oen - 1]][r_[oen]] +
        em().internal_other_mismatch[r_[ien]][r_[ien + 1]][r_[ist - 1]][r_[ist]];

  return energy;
}

Energy Precomp::Hairpin(int st, int en) const {
  const int length = en - st - 1;
  assert(length >= HAIRPIN_MIN_SZ);

  // AU/GU penalty baked into precomputed special table.
  if (length <= MAX_SPECIAL_HAIRPIN_SZ && hairpin[st].special[length] != MAX_E)
    return hairpin[st].special[length];

  const Base stb = r_[st];
  const Base st1b = r_[st + 1];
  const Base en1b = r_[en - 1];
  const Base enb = r_[en];
  Energy energy = em().HairpinInitiation(length) + em().AuGuPenalty(stb, enb);
  const bool all_c = hairpin[st + 1].num_c >= length;

  if (length == 3) {
    if (all_c) energy += em().hairpin_c3_loop;
    return energy;
  }
  energy += em().terminal[r_[st]][st1b][en1b][r_[en]];

  if ((st1b == U && en1b == U) || (st1b == G && en1b == A))
    energy += em().hairpin_uu_ga_first_mismatch;
  if (st1b == G && en1b == G) energy += em().hairpin_gg_first_mismatch;
  if (all_c) energy += em().hairpin_all_c_a * length + em().hairpin_all_c_b;
  if (stb == G && enb == U && st >= 2 && r_[st - 1] == G && r_[st - 2] == G)
    energy += em().hairpin_special_gu_closure;

  return energy;
}

void Precomp::PrecomputeData() {
  // Initialise fast AUGU branch table
  for (Base i = 0; i < 4; ++i)
    for (Base j = 0; j < 4; ++j) augubranch[i][j] = em().multiloop_hack_b + em().AuGuPenalty(i, j);

  const auto min_stack = MinEnergy(&em().stack[0][0][0][0], sizeof(em().stack));

  // Non continuous (-2.1), -4 for WC, -16 for terminal mismatch.
  min_mismatch_coax = em().coax_mismatch_non_contiguous +
      std::min(std::min(em().coax_mismatch_gu_bonus, em().coax_mismatch_wc_bonus), ZERO_E) +
      MinEnergy(&em().terminal[0][0][0][0], sizeof(em().terminal));
  // Minimum of all stacking params.
  min_flush_coax = min_stack;

  Energy min_internal = MinEnergy(&em().internal_1x1[0][0][0][0][0][0], sizeof(em().internal_1x1));
  min_internal = std::min(
      min_internal, MinEnergy(&em().internal_1x2[0][0][0][0][0][0][0], sizeof(em().internal_1x2)));
  min_internal = std::min(min_internal,
      MinEnergy(&em().internal_2x2[0][0][0][0][0][0][0][0], sizeof(em().internal_2x2)));

  verify(em().internal_asym >= ZERO_E,
      "min_internal optimisation does not work for negative asymmetry penalties");
  const auto min_mismatch = 2 *
      std::min(
          MinEnergy(&em().internal_2x3_mismatch[0][0][0][0], sizeof(em().internal_2x3_mismatch)),
          MinEnergy(
              &em().internal_other_mismatch[0][0][0][0], sizeof(em().internal_other_mismatch)));
  const auto min_internal_init = MinEnergy(
      &em().internal_init[4], sizeof(em().internal_init) - 4 * sizeof(em().internal_init[0]));

  // Current implementation needs regular AU/GU for special internal loops.
  const auto min_internal_penalty =
      std::min(2 * std::min(em().internal_au_penalty, em().internal_gu_penalty), ZERO_E);
  const auto min_penalty = std::min(2 * std::min(em().au_penalty, em().gu_penalty), ZERO_E);

  min_internal = std::min(min_internal + min_penalty,
      min_internal_init + min_penalty + min_internal_penalty + min_mismatch);

  const auto min_bulge_init =
      MinEnergy(&em().bulge_init[1], sizeof(em().bulge_init) - sizeof(em().bulge_init[0]));

  const Energy states_bonus = -E(R * T * log(MaxNumContiguous(r_)));
  const Energy min_bulge = min_bulge_init + min_penalty + min_stack +
      std::min(em().bulge_special_c, ZERO_E) + states_bonus;
  min_twoloop_not_stack = std::min(min_bulge, min_internal);

  hairpin = PrecomputeHairpin<HairpinPrecomp<Energy>, /*is_boltz=*/false>(r_, em(), MAX_E);
}

}  // namespace mrna::md::t04
