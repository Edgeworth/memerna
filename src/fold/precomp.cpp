// Copyright 2016, E.
//
// This file is part of memerna.
//
// memerna is free software: you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// memerna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with memerna.
// If not, see <http://www.gnu.org/licenses/>.
#include "precomp.h"
#include "fold/globals.h"
#include "parsing.h"

namespace memerna {
namespace fold {
namespace internal {

using namespace energy;

namespace {

energy_t MinEnergy(const energy_t* energy, std::size_t size) {
  energy_t min = energy[0];
  for (int i = 0; i < int(size / sizeof(energy_t)); ++i)
    min = std::min(min, energy[i]);
  return min;
}
}

int MaxNumContiguous(const primary_t& r) {
  energy_t num_contig = 0;
  energy_t max_num_contig = 0;
  base_t prev = -1;
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

precomp_t PrecomputeData(const primary_t& r, const energy::EnergyModel& em) {
  assert(!r.empty());
  precomp_t pc;
  // Initialise fast AUGU branch table
  for (base_t i = 0; i < 4; ++i)
    for (base_t j = 0; j < 4; ++j)
      pc.augubranch[i][j] = em.multiloop_hack_b + em.AuGuPenalty(i, j);

  const auto min_stack = MinEnergy(&em.stack[0][0][0][0], sizeof(em.stack));

  // Non continuous (-2.1), -4 for WC, -16 for terminal mismatch.
  pc.min_mismatch_coax = em.coax_mismatch_non_contiguous +
      std::min(std::min(em.coax_mismatch_gu_bonus, em.coax_mismatch_wc_bonus), 0) +
      MinEnergy(&em.terminal[0][0][0][0], sizeof(em.terminal));
  // Minimum of all stacking params.
  pc.min_flush_coax = min_stack;

  energy_t min_internal = MinEnergy(&em.internal_1x1[0][0][0][0][0][0], sizeof(em.internal_1x1));
  min_internal = std::min(
      min_internal, MinEnergy(&em.internal_1x2[0][0][0][0][0][0][0], sizeof(em.internal_1x2)));
  min_internal = std::min(
      min_internal, MinEnergy(&em.internal_2x2[0][0][0][0][0][0][0][0], sizeof(em.internal_2x2)));
  verify_expr(em.internal_asym >= 0,
      "min_internal optimisation does not work for negative asymmetry penalties");
  const auto min_mismatch =
      2 *
          std::min(
              MinEnergy(&em.internal_2x3_mismatch[0][0][0][0], sizeof(em.internal_2x3_mismatch)),
              MinEnergy(&em.internal_other_mismatch[0][0][0][0],
                  sizeof(em.internal_other_mismatch)));
  const auto min_internal_init =
      MinEnergy(&em.internal_init[4], sizeof(em.internal_init) - 4 * sizeof(em.internal_init[0]));
  min_internal = std::min(
      min_internal, min_internal_init + std::min(2 * em.internal_augu_penalty, 0) + min_mismatch);

  const auto min_bulge_init =
      MinEnergy(&em.bulge_init[1], sizeof(em.bulge_init) - sizeof(em.bulge_init[0]));

  energy_t states_bonus = -energy_t(round(10.0 * R * T * log(MaxNumContiguous(r))));
  energy_t min_bulge = min_bulge_init + std::min(2 * em.augu_penalty, 0) + min_stack +
      std::min(em.bulge_special_c, 0) + states_bonus;
  pc.min_twoloop_not_stack = std::min(min_bulge, min_internal);

  pc.hairpin.resize(r.size());
  std::string rna_str = parsing::PrimaryToString(r);
  for (const auto& hairpinpair : em.hairpin) {
    const auto& str = hairpinpair.first;
    verify_expr(str.size() - 2 <= hairpin_precomp_t::MAX_SPECIAL_HAIRPIN_SZ,
        "need to increase MAX_SPECIAL_HAIRPIN_SZ");
    auto pos = rna_str.find(str, 0);
    while (pos != std::string::npos) {
      pc.hairpin[pos].special[str.size() - 2] = hairpinpair.second;
      pos = rna_str.find(str, pos + 1);
    }
  }
  const int N = int(r.size());
  pc.hairpin[N - 1].num_c = int(r[N - 1] == C);
  for (int i = N - 2; i >= 0; --i) {
    if (r[i] == C) pc.hairpin[i].num_c = pc.hairpin[i + 1].num_c + 1;
  }

  return pc;
}

energy_t FastTwoLoop(int ost, int oen, int ist, int ien) {
  int toplen = ist - ost - 1, botlen = oen - ien - 1;
  if (toplen == 0 && botlen == 0) return gem.stack[gr[ost]][gr[ist]][gr[ien]][gr[oen]];
  if (toplen == 0 || botlen == 0) return gem.Bulge(gr, ost, oen, ist, ien);
  if (toplen == 1 && botlen == 1)
    return gem.internal_1x1[gr[ost]][gr[ost + 1]][gr[ist]][gr[ien]][gr[ien + 1]][gr[oen]];
  if (toplen == 1 && botlen == 2)
    return gem.internal_1x2[gr[ost]][gr[ost + 1]][gr[ist]][gr[ien]][gr[ien + 1]][gr[ien + 2]]
    [gr[oen]];
  if (toplen == 2 && botlen == 1)
    return gem.internal_1x2[gr[ien]][gr[ien + 1]][gr[oen]][gr[ost]][gr[ost + 1]][gr[ost + 2]]
    [gr[ist]];
  if (toplen == 2 && botlen == 2)
    return gem.internal_2x2[gr[ost]][gr[ost + 1]][gr[ost + 2]][gr[ist]][gr[ien]][gr[ien + 1]]
    [gr[ien + 2]][gr[oen]];

  static_assert(
      TWOLOOP_MAX_SZ <= EnergyModel::INITIATION_CACHE_SZ, "initiation cache not large enough");
  energy_t energy = gem.internal_init[toplen + botlen] +
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

energy_t FastHairpin(int st, int en) {
  int length = en - st - 1;
  assert(length >= HAIRPIN_MIN_SZ);
  if (length <= internal::hairpin_precomp_t::MAX_SPECIAL_HAIRPIN_SZ &&
      gpc.hairpin[st].special[length] != MAX_E)
    return gpc.hairpin[st].special[length];
  base_t stb = gr[st], st1b = gr[st + 1], en1b = gr[en - 1], enb = gr[en];
  energy_t energy = gem.HairpinInitiation(length) + gem.AuGuPenalty(stb, enb);

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
}
}
}
