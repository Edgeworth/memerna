// Copyright 2016 Eliot Courtney.
#include "compute/energy/fast_energy.h"

#include <algorithm>
#include <cassert>

#include "compute/energy/globals.h"
#include "model/parsing.h"
#include "util/macros.h"

namespace mrna::energy {

namespace {

Energy MinEnergy(const Energy* energy, std::size_t size) {
  Energy min = energy[0];
  for (int i = 0; i < static_cast<int>(size / sizeof(Energy)); ++i) min = std::min(min, energy[i]);
  return min;
}

}  // namespace

namespace internal {

int MaxNumContiguous(const Primary& r) {
  Energy num_contig = 0;
  Energy max_num_contig = 0;
  Base prev = -1;
  for (auto b : r) {
    if (b == prev)
      num_contig++;
    else
      num_contig = 1;
    prev = b;
    max_num_contig = std::max(max_num_contig, num_contig);
  }
  return max_num_contig;
}

}  // namespace internal

Precomp PrecomputeData(const Primary& r, const energy::EnergyModel& em) {
  assert(!r.empty());
  Precomp pc;
  // Initialise fast AUGU branch table
  for (Base i = 0; i < 4; ++i)
    for (Base j = 0; j < 4; ++j) pc.augubranch[i][j] = em.multiloop_hack_b + em.AuGuPenalty(i, j);

  const auto min_stack = MinEnergy(&em.stack[0][0][0][0], sizeof(em.stack));

  // Non continuous (-2.1), -4 for WC, -16 for terminal mismatch.
  pc.min_mismatch_coax = em.coax_mismatch_non_contiguous +
      std::min(std::min(em.coax_mismatch_gu_bonus, em.coax_mismatch_wc_bonus), 0) +
      MinEnergy(&em.terminal[0][0][0][0], sizeof(em.terminal));
  // Minimum of all stacking params.
  pc.min_flush_coax = min_stack;

  Energy min_internal = MinEnergy(&em.internal_1x1[0][0][0][0][0][0], sizeof(em.internal_1x1));
  min_internal = std::min(
      min_internal, MinEnergy(&em.internal_1x2[0][0][0][0][0][0][0], sizeof(em.internal_1x2)));
  min_internal = std::min(
      min_internal, MinEnergy(&em.internal_2x2[0][0][0][0][0][0][0][0], sizeof(em.internal_2x2)));
  verify(em.internal_asym >= 0,
      "min_internal optimisation does not work for negative asymmetry penalties");
  const auto min_mismatch = 2 *
      std::min(MinEnergy(&em.internal_2x3_mismatch[0][0][0][0], sizeof(em.internal_2x3_mismatch)),
          MinEnergy(&em.internal_other_mismatch[0][0][0][0], sizeof(em.internal_other_mismatch)));
  const auto min_internal_init =
      MinEnergy(&em.internal_init[4], sizeof(em.internal_init) - 4 * sizeof(em.internal_init[0]));
  min_internal = std::min(
      min_internal, min_internal_init + std::min(2 * em.internal_augu_penalty, 0) + min_mismatch);

  const auto min_bulge_init =
      MinEnergy(&em.bulge_init[1], sizeof(em.bulge_init) - sizeof(em.bulge_init[0]));

  Energy states_bonus = -Energy(round(10.0 * R * T * log(internal::MaxNumContiguous(r))));
  Energy min_bulge = min_bulge_init + std::min(2 * em.augu_penalty, 0) + min_stack +
      std::min(em.bulge_special_c, 0) + states_bonus;
  pc.min_twoloop_not_stack = std::min(min_bulge, min_internal);

  pc.hairpin = PrecomputeHairpin<HairpinPrecomp<Energy, MAX_E>>(r, em);

  return pc;
}

Energy FastTwoLoop(int ost, int oen, int ist, int ien) {
  int toplen = ist - ost - 1, botlen = oen - ien - 1;
  if (toplen == 0 && botlen == 0) return gem.stack[gr[ost]][gr[ist]][gr[ien]][gr[oen]];
  if (toplen == 0 || botlen == 0) return gem.Bulge(gr, ost, oen, ist, ien);
  if (toplen == 1 && botlen == 1)
    return gem.internal_1x1[gr[ost]][gr[ost + 1]][gr[ist]][gr[ien]][gr[ien + 1]][gr[oen]];
  if (toplen == 1 && botlen == 2)
    return gem
        .internal_1x2[gr[ost]][gr[ost + 1]][gr[ist]][gr[ien]][gr[ien + 1]][gr[ien + 2]][gr[oen]];
  if (toplen == 2 && botlen == 1)
    return gem
        .internal_1x2[gr[ien]][gr[ien + 1]][gr[oen]][gr[ost]][gr[ost + 1]][gr[ost + 2]][gr[ist]];
  if (toplen == 2 && botlen == 2)
    return gem.internal_2x2[gr[ost]][gr[ost + 1]][gr[ost + 2]][gr[ist]][gr[ien]][gr[ien + 1]]
                           [gr[ien + 2]][gr[oen]];

  static_assert(
      TWOLOOP_MAX_SZ <= EnergyModel::INITIATION_CACHE_SZ, "initiation cache not large enough");
  assert(toplen + botlen < EnergyModel::INITIATION_CACHE_SZ);
  Energy energy = gem.internal_init[toplen + botlen] +
      std::min(std::abs(toplen - botlen) * gem.internal_asym, NINIO_MAX_ASYM);

  energy += gem.InternalLoopAuGuPenalty(gr[ost], gr[oen]);
  energy += gem.InternalLoopAuGuPenalty(gr[ist], gr[ien]);

  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2))
    energy += gem.internal_2x3_mismatch[gr[ost]][gr[ost + 1]][gr[oen - 1]][gr[oen]] +
        gem.internal_2x3_mismatch[gr[ien]][gr[ien + 1]][gr[ist - 1]][gr[ist]];
  else if (toplen != 1 && botlen != 1)
    energy += gem.internal_other_mismatch[gr[ost]][gr[ost + 1]][gr[oen - 1]][gr[oen]] +
        gem.internal_other_mismatch[gr[ien]][gr[ien + 1]][gr[ist - 1]][gr[ist]];

  return energy;
}

Energy FastHairpin(int st, int en) {
  int length = en - st - 1;
  assert(length >= HAIRPIN_MIN_SZ);
  if (length <= MAX_SPECIAL_HAIRPIN_SZ && gpc.hairpin[st].special[length] != MAX_E)
    return gpc.hairpin[st].special[length];
  Base stb = gr[st], st1b = gr[st + 1], en1b = gr[en - 1], enb = gr[en];
  Energy energy = gem.HairpinInitiation(length) + gem.AuGuPenalty(stb, enb);

  bool all_c = gpc.hairpin[st + 1].num_c >= length;

  if (length == 3) {
    if (all_c) energy += gem.hairpin_c3_loop;
    return energy;
  }
  energy += gem.terminal[gr[st]][st1b][en1b][gr[en]];

  if ((st1b == U && en1b == U) || (st1b == G && en1b == A))
    energy += gem.hairpin_uu_ga_first_mismatch;
  if (st1b == G && en1b == G) energy += gem.hairpin_gg_first_mismatch;
  if (all_c) energy += gem.hairpin_all_c_a * length + gem.hairpin_all_c_b;
  if (stb == G && enb == U && st >= 2 && gr[st - 1] == G && gr[st - 2] == G)
    energy += gem.hairpin_special_gu_closure;

  return energy;
}

}  // namespace mrna::energy
