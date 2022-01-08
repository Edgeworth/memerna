// Copyright 2016 Eliot Courtney.
#include "compute/energy/precomp.h"

#include <algorithm>
#include <cassert>
#include <memory>
#include <string>

#include "compute/energy/structure.h"
#include "model/parsing.h"
#include "util/string.h"

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

Precomp::Precomp(Primary r, EnergyModel em) : r_(r), em_(em) { PrecomputeData(); }

Energy Precomp::FastTwoLoop(int ost, int oen, int ist, int ien) const {
  int toplen = ist - ost - 1, botlen = oen - ien - 1;
  if (toplen == 0 && botlen == 0) return em_.stack[r_[ost]][r_[ist]][r_[ien]][r_[oen]];
  if (toplen == 0 || botlen == 0) return em_.Bulge(r_, ost, oen, ist, ien);
  if (toplen == 1 && botlen == 1)
    return em_.internal_1x1[r_[ost]][r_[ost + 1]][r_[ist]][r_[ien]][r_[ien + 1]][r_[oen]];
  if (toplen == 1 && botlen == 2)
    return em_
        .internal_1x2[r_[ost]][r_[ost + 1]][r_[ist]][r_[ien]][r_[ien + 1]][r_[ien + 2]][r_[oen]];
  if (toplen == 2 && botlen == 1)
    return em_
        .internal_1x2[r_[ien]][r_[ien + 1]][r_[oen]][r_[ost]][r_[ost + 1]][r_[ost + 2]][r_[ist]];
  if (toplen == 2 && botlen == 2)
    return em_.internal_2x2[r_[ost]][r_[ost + 1]][r_[ost + 2]][r_[ist]][r_[ien]][r_[ien + 1]]
                           [r_[ien + 2]][r_[oen]];

  static_assert(
      TWOLOOP_MAX_SZ <= EnergyModel::INITIATION_CACHE_SZ, "initiation cache not large enough");
  assert(toplen + botlen < EnergyModel::INITIATION_CACHE_SZ);
  Energy energy = em_.internal_init[toplen + botlen] +
      std::min(std::abs(toplen - botlen) * em_.internal_asym, NINIO_MAX_ASYM);

  energy += em_.InternalLoopAuGuPenalty(r_[ost], r_[oen]);
  energy += em_.InternalLoopAuGuPenalty(r_[ist], r_[ien]);

  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2))
    energy += em_.internal_2x3_mismatch[r_[ost]][r_[ost + 1]][r_[oen - 1]][r_[oen]] +
        em_.internal_2x3_mismatch[r_[ien]][r_[ien + 1]][r_[ist - 1]][r_[ist]];
  else if (toplen != 1 && botlen != 1)
    energy += em_.internal_other_mismatch[r_[ost]][r_[ost + 1]][r_[oen - 1]][r_[oen]] +
        em_.internal_other_mismatch[r_[ien]][r_[ien + 1]][r_[ist - 1]][r_[ist]];

  return energy;
}

Energy Precomp::FastHairpin(int st, int en) const {
  int length = en - st - 1;
  assert(length >= HAIRPIN_MIN_SZ);
  if (length <= MAX_SPECIAL_HAIRPIN_SZ && hairpin[st].special[length] != MAX_E)
    return hairpin[st].special[length];
  Base stb = r_[st], st1b = r_[st + 1], en1b = r_[en - 1], enb = r_[en];
  Energy energy = em_.HairpinInitiation(length) + em_.AuGuPenalty(stb, enb);

  bool all_c = hairpin[st + 1].num_c >= length;

  if (length == 3) {
    if (all_c) energy += em_.hairpin_c3_loop;
    return energy;
  }
  energy += em_.terminal[r_[st]][st1b][en1b][r_[en]];

  if ((st1b == U && en1b == U) || (st1b == G && en1b == A))
    energy += em_.hairpin_uu_ga_first_mismatch;
  if (st1b == G && en1b == G) energy += em_.hairpin_gg_first_mismatch;
  if (all_c) energy += em_.hairpin_all_c_a * length + em_.hairpin_all_c_b;
  if (stb == G && enb == U && st >= 2 && r_[st - 1] == G && r_[st - 2] == G)
    energy += em_.hairpin_special_gu_closure;

  return energy;
}

void Precomp::PrecomputeData() {
  // Initialise fast AUGU branch table
  for (Base i = 0; i < 4; ++i)
    for (Base j = 0; j < 4; ++j) augubranch[i][j] = em_.multiloop_hack_b + em_.AuGuPenalty(i, j);

  const auto min_stack = MinEnergy(&em_.stack[0][0][0][0], sizeof(em_.stack));

  // Non continuous (-2.1), -4 for WC, -16 for terminal mismatch.
  min_mismatch_coax = em_.coax_mismatch_non_contiguous +
      std::min(std::min(em_.coax_mismatch_gu_bonus, em_.coax_mismatch_wc_bonus), 0) +
      MinEnergy(&em_.terminal[0][0][0][0], sizeof(em_.terminal));
  // Minimum of all stacking params.
  min_flush_coax = min_stack;

  Energy min_internal = MinEnergy(&em_.internal_1x1[0][0][0][0][0][0], sizeof(em_.internal_1x1));
  min_internal = std::min(
      min_internal, MinEnergy(&em_.internal_1x2[0][0][0][0][0][0][0], sizeof(em_.internal_1x2)));
  min_internal = std::min(
      min_internal, MinEnergy(&em_.internal_2x2[0][0][0][0][0][0][0][0], sizeof(em_.internal_2x2)));
  verify(em_.internal_asym >= 0,
      "min_internal optimisation does not work for negative asymmetry penalties");
  const auto min_mismatch = 2 *
      std::min(MinEnergy(&em_.internal_2x3_mismatch[0][0][0][0], sizeof(em_.internal_2x3_mismatch)),
          MinEnergy(&em_.internal_other_mismatch[0][0][0][0], sizeof(em_.internal_other_mismatch)));
  const auto min_internal_init = MinEnergy(
      &em_.internal_init[4], sizeof(em_.internal_init) - 4 * sizeof(em_.internal_init[0]));
  min_internal = std::min(
      min_internal, min_internal_init + std::min(2 * em_.internal_augu_penalty, 0) + min_mismatch);

  const auto min_bulge_init =
      MinEnergy(&em_.bulge_init[1], sizeof(em_.bulge_init) - sizeof(em_.bulge_init[0]));

  Energy states_bonus = -Energy(round(10.0 * R * T * log(internal::MaxNumContiguous(r_))));
  Energy min_bulge = min_bulge_init + std::min(2 * em_.augu_penalty, 0) + min_stack +
      std::min(em_.bulge_special_c, 0) + states_bonus;
  min_twoloop_not_stack = std::min(min_bulge, min_internal);

  hairpin = PrecomputeHairpin<HairpinPrecomp<Energy, MAX_E>>(r_, em_);
}

}  // namespace mrna::energy
