// Copyright 2022 Eliot Courtney.
#include "backends/baseopt/energy/precomp.h"

namespace mrna::md::base::opt {

Precomp::Precomp(Primary r, Model::Ptr m) : PrecompBase(std::move(r), std::move(m)) {}

Energy Precomp::TwoLoop(int ost, int oen, int ist, int ien) const {
  const int toplen = ist - ost - 1;
  const int botlen = oen - ien - 1;
  if (toplen == 0 && botlen == 0) return m_->stack[r_[ost]][r_[ist]][r_[ien]][r_[oen]];
  if (toplen == 0 || botlen == 0) return m_->Bulge(r_, ost, oen, ist, ien);

  Energy energy = ZERO_E;
  energy += m_->AuGuPenalty(r_[ost], r_[oen]);
  energy += m_->AuGuPenalty(r_[ist], r_[ien]);

  if (toplen == 1 && botlen == 1)
    return energy + m_->internal_1x1[r_[ost]][r_[ost + 1]][r_[ist]][r_[ien]][r_[ien + 1]][r_[oen]];
  if (toplen == 1 && botlen == 2)
    return energy +
        m_->internal_1x2[r_[ost]][r_[ost + 1]][r_[ist]][r_[ien]][r_[ien + 1]][r_[ien + 2]][r_[oen]];
  if (toplen == 2 && botlen == 1)
    return energy +
        m_->internal_1x2[r_[ien]][r_[ien + 1]][r_[oen]][r_[ost]][r_[ost + 1]][r_[ost + 2]][r_[ist]];
  if (toplen == 2 && botlen == 2)
    return energy +
        m_->internal_2x2[r_[ost]][r_[ost + 1]][r_[ost + 2]][r_[ist]][r_[ien]][r_[ien + 1]]
                        [r_[ien + 2]][r_[oen]];

  static_assert(
      TWOLOOP_MAX_SZ <= ModelBase::INITIATION_CACHE_SZ, "initiation cache not large enough");
  assert(toplen + botlen < ModelBase::INITIATION_CACHE_SZ);
  energy += m_->internal_init[toplen + botlen] +
      std::min(std::abs(toplen - botlen) * m_->internal_asym, NINIO_MAX_ASYM);

  energy += m_->InternalLoopAuGuPenalty(r_[ost], r_[oen]);
  energy += m_->InternalLoopAuGuPenalty(r_[ist], r_[ien]);

  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2))
    energy += m_->internal_2x3_mismatch[r_[ost]][r_[ost + 1]][r_[oen - 1]][r_[oen]] +
        m_->internal_2x3_mismatch[r_[ien]][r_[ien + 1]][r_[ist - 1]][r_[ist]];
  else if (toplen != 1 && botlen != 1)
    energy += m_->internal_other_mismatch[r_[ost]][r_[ost + 1]][r_[oen - 1]][r_[oen]] +
        m_->internal_other_mismatch[r_[ien]][r_[ien + 1]][r_[ist - 1]][r_[ist]];

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
  Energy energy = m_->HairpinInitiation(length) + m_->AuGuPenalty(stb, enb);
  const bool all_c = hairpin[st + 1].num_c >= length;

  if (length == 3) {
    if (all_c) energy += m_->hairpin_c3_loop;
    return energy;
  }
  energy += m_->terminal[r_[st]][st1b][en1b][r_[en]];

  if ((st1b == U && en1b == U) || (st1b == G && en1b == A))
    energy += m_->hairpin_uu_ga_first_mismatch;
  if (st1b == G && en1b == G) energy += m_->hairpin_gg_first_mismatch;
  if (all_c) energy += m_->hairpin_all_c_a * length + m_->hairpin_all_c_b;
  if (stb == G && enb == U && st >= 2 && r_[st - 1] == G && r_[st - 2] == G)
    energy += m_->hairpin_special_gu_closure;

  return energy;
}

}  // namespace mrna::md::base::opt
