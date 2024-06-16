// Copyright 2022 Eliot Courtney.
#include "backends/baseopt/energy/boltz_precomp.h"

namespace mrna::md::base::opt {

BoltzPrecomp::BoltzPrecomp(Primary r, BoltzModel::Ptr bm)
    : BoltzPrecompBase(std::move(r), std::move(bm)) {}

BoltzEnergy BoltzPrecomp::Hairpin(int st, int en) const {
  const auto& m = bm_->m();
  const int length = en - st - 1;
  assert(length >= HAIRPIN_MIN_SZ);
  if (length <= MAX_SPECIAL_HAIRPIN_SZ && hairpin[st].special[length] > -1)
    return hairpin[st].special[length];
  const Base stb = r_[st];
  const Base st1b = r_[st + 1];
  const Base en1b = r_[en - 1];
  const Base enb = r_[en];

  BoltzEnergy energy = (m.HairpinInitiation(length) + m.AuGuPenalty(stb, enb)).Boltz();

  const bool all_c = hairpin[st + 1].num_c >= length;

  if (length == 3) {
    if (all_c) energy *= bm_->hairpin_c3_loop;
    return energy;
  }
  energy *= bm_->terminal[r_[st]][st1b][en1b][r_[en]];

  if ((st1b == U && en1b == U) || (st1b == G && en1b == A))
    energy *= bm_->hairpin_uu_ga_first_mismatch;
  if (st1b == G && en1b == G) energy *= bm_->hairpin_gg_first_mismatch;
  if (all_c) energy *= (m.hairpin_all_c_a * length).Boltz() * bm_->hairpin_all_c_b;
  if (stb == G && enb == U && st >= 2 && r_[st - 1] == G && r_[st - 2] == G)
    energy *= bm_->hairpin_special_gu_closure;

  return energy;
}

BoltzEnergy BoltzPrecomp::TwoLoop(int ost, int oen, int ist, int ien) const {
  const auto& m = bm_->m();
  const int toplen = ist - ost - 1;
  const int botlen = oen - ien - 1;
  if (toplen == 0 && botlen == 0) return bm_->stack[r_[ost]][r_[ist]][r_[ien]][r_[oen]];
  if (toplen == 0 || botlen == 0) return bm_->Bulge(r_, ost, oen, ist, ien);

  BoltzEnergy energy = bm_->AuGuPenalty(r_[ost], r_[oen]) * bm_->AuGuPenalty(r_[ist], r_[ien]);

  if (toplen == 1 && botlen == 1)
    return energy * bm_->internal_1x1[r_[ost]][r_[ost + 1]][r_[ist]][r_[ien]][r_[ien + 1]][r_[oen]];
  if (toplen == 1 && botlen == 2)
    return energy *
        bm_->internal_1x2[r_[ost]][r_[ost + 1]][r_[ist]][r_[ien]][r_[ien + 1]][r_[ien + 2]]
                         [r_[oen]];
  if (toplen == 2 && botlen == 1)
    return energy *
        bm_->internal_1x2[r_[ien]][r_[ien + 1]][r_[oen]][r_[ost]][r_[ost + 1]][r_[ost + 2]]
                         [r_[ist]];
  if (toplen == 2 && botlen == 2)
    return energy *
        bm_->internal_2x2[r_[ost]][r_[ost + 1]][r_[ost + 2]][r_[ist]][r_[ien]][r_[ien + 1]]
                         [r_[ien + 2]][r_[oen]];

  static_assert(TWOLOOP_MAX_SZ <= std::remove_reference_t<decltype(m)>::INITIATION_CACHE_SZ,
      "initiation cache not large enough");
  assert(toplen + botlen < std::remove_reference_t<decltype(m)>::INITIATION_CACHE_SZ);
  energy *= bm_->internal_init[toplen + botlen] *
      std::min(std::abs(toplen - botlen) * m.internal_asym, NINIO_MAX_ASYM).Boltz();

  energy *= bm_->InternalLoopAuGuPenalty(r_[ost], r_[oen]);
  energy *= bm_->InternalLoopAuGuPenalty(r_[ist], r_[ien]);

  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2))
    energy *= bm_->internal_2x3_mismatch[r_[ost]][r_[ost + 1]][r_[oen - 1]][r_[oen]] *
        bm_->internal_2x3_mismatch[r_[ien]][r_[ien + 1]][r_[ist - 1]][r_[ist]];
  else if (toplen != 1 && botlen != 1)
    energy *= bm_->internal_other_mismatch[r_[ost]][r_[ost + 1]][r_[oen - 1]][r_[oen]] *
        bm_->internal_other_mismatch[r_[ien]][r_[ien + 1]][r_[ist - 1]][r_[ist]];

  return energy;
}

}  // namespace mrna::md::base::opt
