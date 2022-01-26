// Copyright 2022 E.
#include "compute/energy/boltzmann_precomp.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <memory>
#include <unordered_map>
#include <utility>

#include "compute/energy/model.h"
#include "model/base.h"

namespace mrna::energy {

BoltzPrecomp::BoltzPrecomp(Primary r, BoltzEnergyModel bem)
    : r_(std::move(r)), bem_(BoltzEnergyModel(std::move(bem))) {
  PrecomputeData();
}

BoltzEnergy BoltzPrecomp::Hairpin(int st, int en) const {
  int length = en - st - 1;
  assert(length >= HAIRPIN_MIN_SZ);
  if (length <= MAX_SPECIAL_HAIRPIN_SZ && hairpin[st].special[length] > -1)
    return hairpin[st].special[length];
  Base stb = r_[st], st1b = r_[st + 1], en1b = r_[en - 1], enb = r_[en];
  BoltzEnergy energy = Boltz(bem_.em().HairpinInitiation(length) + bem_.em().AuGuPenalty(stb, enb));

  const bool all_c = hairpin[st + 1].num_c >= length;

  if (length == 3) {
    if (all_c) energy *= bem_.hairpin_c3_loop;
    return energy;
  }
  energy *= bem_.terminal[r_[st]][st1b][en1b][r_[en]];

  if ((st1b == U && en1b == U) || (st1b == G && en1b == A))
    energy *= bem_.hairpin_uu_ga_first_mismatch;
  if (st1b == G && en1b == G) energy *= bem_.hairpin_gg_first_mismatch;
  if (all_c) energy *= Boltz(bem_.em().hairpin_all_c_a * length) * bem_.hairpin_all_c_b;
  if (stb == G && enb == U && st >= 2 && r_[st - 1] == G && r_[st - 2] == G)
    energy *= bem_.hairpin_special_gu_closure;

  return energy;
}

BoltzEnergy BoltzPrecomp::TwoLoop(int ost, int oen, int ist, int ien) const {
  const int toplen = ist - ost - 1, botlen = oen - ien - 1;
  if (toplen == 0 && botlen == 0) return bem_.stack[r_[ost]][r_[ist]][r_[ien]][r_[oen]];
  if (toplen == 0 || botlen == 0) return bem_.Bulge(r_, ost, oen, ist, ien);
  if (toplen == 1 && botlen == 1)
    return bem_.internal_1x1[r_[ost]][r_[ost + 1]][r_[ist]][r_[ien]][r_[ien + 1]][r_[oen]];
  if (toplen == 1 && botlen == 2)
    return bem_
        .internal_1x2[r_[ost]][r_[ost + 1]][r_[ist]][r_[ien]][r_[ien + 1]][r_[ien + 2]][r_[oen]];
  if (toplen == 2 && botlen == 1)
    return bem_
        .internal_1x2[r_[ien]][r_[ien + 1]][r_[oen]][r_[ost]][r_[ost + 1]][r_[ost + 2]][r_[ist]];
  if (toplen == 2 && botlen == 2)
    return bem_.internal_2x2[r_[ost]][r_[ost + 1]][r_[ost + 2]][r_[ist]][r_[ien]][r_[ien + 1]]
                            [r_[ien + 2]][r_[oen]];

  static_assert(
      TWOLOOP_MAX_SZ <= EnergyModel::INITIATION_CACHE_SZ, "initiation cache not large enough");
  assert(toplen + botlen < EnergyModel::INITIATION_CACHE_SZ);
  BoltzEnergy energy = bem_.internal_init[toplen + botlen] *
      Boltz(std::min(std::abs(toplen - botlen) * bem_.em().internal_asym, NINIO_MAX_ASYM));

  energy *= bem_.InternalLoopAuGuPenalty(r_[ost], r_[oen]);
  energy *= bem_.InternalLoopAuGuPenalty(r_[ist], r_[ien]);

  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2))
    energy *= bem_.internal_2x3_mismatch[r_[ost]][r_[ost + 1]][r_[oen - 1]][r_[oen]] *
        bem_.internal_2x3_mismatch[r_[ien]][r_[ien + 1]][r_[ist - 1]][r_[ist]];
  else if (toplen != 1 && botlen != 1)
    energy *= bem_.internal_other_mismatch[r_[ost]][r_[ost + 1]][r_[oen - 1]][r_[oen]] *
        bem_.internal_other_mismatch[r_[ien]][r_[ien + 1]][r_[ist - 1]][r_[ist]];

  return energy;
}

void BoltzPrecomp::PrecomputeData() {
  for (Base i = 0; i < 4; ++i)
    for (Base j = 0; j < 4; ++j) augubranch[i][j] = bem_.multiloop_hack_b * bem_.AuGuPenalty(i, j);
  hairpin = PrecomputeHairpin<HairpinPrecomp<BoltzEnergy, -1>>(r_, bem_);
}

}  // namespace mrna::energy