// Copyright 2022 Eliot Courtney.
#include "compute/energy/t04/boltz_precomp.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <memory>
#include <unordered_map>
#include <utility>

#include "compute/energy/energy.h"
#include "model/base.h"

namespace mrna::energy::t04 {

BoltzPrecomp::BoltzPrecomp(Primary r, BoltzModelPtr bem) : r_(std::move(r)), bem_(std::move(bem)) {
  PrecomputeData();
}

BoltzEnergy BoltzPrecomp::Hairpin(int st, int en) const {
  int length = en - st - 1;
  assert(length >= HAIRPIN_MIN_SZ);
  if (length <= MAX_SPECIAL_HAIRPIN_SZ && hairpin[st].special[length] > -1)
    return hairpin[st].special[length];
  Base stb = r_[st];
  Base st1b = r_[st + 1];
  Base en1b = r_[en - 1];
  Base enb = r_[en];
  BoltzEnergy energy = (em()->HairpinInitiation(length) + em()->AuGuPenalty(stb, enb)).Boltz();

  const bool all_c = hairpin[st + 1].num_c >= length;

  if (length == 3) {
    if (all_c) energy *= bem().hairpin_c3_loop;
    return energy;
  }
  energy *= bem().terminal[r_[st]][st1b][en1b][r_[en]];

  if ((st1b == U && en1b == U) || (st1b == G && en1b == A))
    energy *= bem().hairpin_uu_ga_first_mismatch;
  if (st1b == G && en1b == G) energy *= bem().hairpin_gg_first_mismatch;
  if (all_c) energy *= (em()->hairpin_all_c_a * length).Boltz() * bem().hairpin_all_c_b;
  if (stb == G && enb == U && st >= 2 && r_[st - 1] == G && r_[st - 2] == G)
    energy *= bem().hairpin_special_gu_closure;

  return energy;
}

BoltzEnergy BoltzPrecomp::TwoLoop(int ost, int oen, int ist, int ien) const {
  const int toplen = ist - ost - 1;
  const int botlen = oen - ien - 1;
  if (toplen == 0 && botlen == 0) return bem().stack[r_[ost]][r_[ist]][r_[ien]][r_[oen]];
  if (toplen == 0 || botlen == 0) return bem().Bulge(r_, ost, oen, ist, ien);
  if (toplen == 1 && botlen == 1)
    return bem().internal_1x1[r_[ost]][r_[ost + 1]][r_[ist]][r_[ien]][r_[ien + 1]][r_[oen]];
  if (toplen == 1 && botlen == 2)
    return bem()
        .internal_1x2[r_[ost]][r_[ost + 1]][r_[ist]][r_[ien]][r_[ien + 1]][r_[ien + 2]][r_[oen]];
  if (toplen == 2 && botlen == 1)
    return bem()
        .internal_1x2[r_[ien]][r_[ien + 1]][r_[oen]][r_[ost]][r_[ost + 1]][r_[ost + 2]][r_[ist]];
  if (toplen == 2 && botlen == 2)
    return bem().internal_2x2[r_[ost]][r_[ost + 1]][r_[ost + 2]][r_[ist]][r_[ien]][r_[ien + 1]]
                             [r_[ien + 2]][r_[oen]];

  static_assert(TWOLOOP_MAX_SZ <= Model::INITIATION_CACHE_SZ, "initiation cache not large enough");
  assert(toplen + botlen < Model::INITIATION_CACHE_SZ);
  BoltzEnergy energy = bem().internal_init[toplen + botlen] *
      std::min(std::abs(toplen - botlen) * em()->internal_asym, NINIO_MAX_ASYM).Boltz();

  energy *= bem().InternalLoopAuGuPenalty(r_[ost], r_[oen]);
  energy *= bem().InternalLoopAuGuPenalty(r_[ist], r_[ien]);

  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2))
    energy *= bem().internal_2x3_mismatch[r_[ost]][r_[ost + 1]][r_[oen - 1]][r_[oen]] *
        bem().internal_2x3_mismatch[r_[ien]][r_[ien + 1]][r_[ist - 1]][r_[ist]];
  else if (toplen != 1 && botlen != 1)
    energy *= bem().internal_other_mismatch[r_[ost]][r_[ost + 1]][r_[oen - 1]][r_[oen]] *
        bem().internal_other_mismatch[r_[ien]][r_[ien + 1]][r_[ist - 1]][r_[ist]];

  return energy;
}

void BoltzPrecomp::PrecomputeData() {
  for (Base i = 0; i < 4; ++i)
    for (Base j = 0; j < 4; ++j)
      augubranch[i][j] = bem().multiloop_hack_b * bem().AuGuPenalty(i, j);
  hairpin = PrecomputeHairpin<HairpinPrecomp<BoltzEnergy>>(r_, bem(), -1.0);
}

}  // namespace mrna::energy::t04
