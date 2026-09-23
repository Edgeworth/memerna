// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_BASE_ENERGY_PRECOMP_H_
#define BACKENDS_BASE_ENERGY_PRECOMP_H_

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <utility>

#include "api/energy/energy_cfg.h"
#include "backends/base/energy/model.h"
#include "backends/common/base/precomp_base.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::base {

template <erg::EnergyCfg CFG, class PseudofreeCfgT, class ModelT = Model>
class Precomp : public PrecompBase<ModelT> {
 public:
  PseudofreeCfgT pf;
  Energy min_pf_unpaired{};
  Energy sum_neg_pf{};

  Precomp(Primary r, ModelT::Ptr m, PseudofreeCfgT pf);

  [[nodiscard]] Energy TwoLoop(int ost, int oen, int ist, int ien) const;

  [[nodiscard]] Energy Hairpin(int st, int en) const;
};

template <erg::EnergyCfg CFG, class PseudofreeCfgT, class ModelT>
Precomp<CFG, PseudofreeCfgT, ModelT>::Precomp(Primary r, ModelT::Ptr m, PseudofreeCfgT pf_)
    : PrecompBase<ModelT>(std::move(r), std::move(m)), pf(std::move(pf_)),
      min_pf_unpaired(pf.MinUnpaired()), sum_neg_pf(pf.SmallestSubsetSum()) {
  pf.Verify(this->r_);
}

template <erg::EnergyCfg CFG, class PseudofreeCfgT, class ModelT>
Energy Precomp<CFG, PseudofreeCfgT, ModelT>::TwoLoop(int ost, int oen, int ist, int ien) const {
  const auto& r = this->r_;
  const auto& m = this->m_;
  const int toplen = ist - ost - 1;
  const int botlen = oen - ien - 1;

  if (toplen == 0 && botlen == 0)
    return pf.Paired(ost, oen) + m->stack[r[ost]][r[ist]][r[ien]][r[oen]];
  if (toplen == 0 || botlen == 0) return m->template Bulge<CFG>(r, pf, ost, oen, ist, ien);

  Energy energy = m->AuGuPenalty(r[ost], r[oen]) + m->AuGuPenalty(r[ist], r[ien]) +
      pf.Paired(ost, oen) + pf.UnpairedSum(ost + 1, ist - 1) + pf.UnpairedSum(ien + 1, oen - 1);

  if (toplen == 1 && botlen == 1)
    return energy + m->internal_1x1[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[oen]];
  if (toplen == 1 && botlen == 2)
    return energy +
        m->internal_1x2[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]][r[oen]];
  if (toplen == 2 && botlen == 1)
    return energy +
        m->internal_1x2[r[ien]][r[ien + 1]][r[oen]][r[ost]][r[ost + 1]][r[ost + 2]][r[ist]];
  if (toplen == 2 && botlen == 2)
    return energy +
        m->internal_2x2[r[ost]][r[ost + 1]][r[ost + 2]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]]
                       [r[oen]];

  static_assert(TWOLOOP_MAX_SZ <= ModelT::INITIATION_CACHE_SZ, "initiation cache not large enough");
  assert(toplen + botlen < ModelT::INITIATION_CACHE_SZ);
  energy += m->internal_init[toplen + botlen] +
      std::min(std::abs(toplen - botlen) * m->internal_asym, NINIO_MAX_ASYM);

  energy += m->InternalLoopAuGuPenalty(r[ost], r[oen]);
  energy += m->InternalLoopAuGuPenalty(r[ist], r[ien]);

  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2))
    energy += m->internal_2x3_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
        m->internal_2x3_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];
  else if (toplen != 1 && botlen != 1)
    energy += m->internal_other_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
        m->internal_other_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];

  return energy;
}

template <erg::EnergyCfg CFG, class PseudofreeCfgT, class ModelT>
Energy Precomp<CFG, PseudofreeCfgT, ModelT>::Hairpin(int st, int en) const {
  const auto& r = this->r_;
  const auto& m = this->m_;
  const auto& hairpin = this->hairpin;
  const int length = en - st - 1;
  assert(length >= HAIRPIN_MIN_SZ);

  Energy energy = pf.UnpairedSum(st + 1, en - 1) + pf.Paired(st, en);

  // AU/GU penalty baked into precomputed special table.
  if (length <= MAX_SPECIAL_HAIRPIN_SZ && hairpin[st].special[length] != MAX_E)
    return energy + hairpin[st].special[length];

  const Base stb = r[st];
  const Base st1b = r[st + 1];
  const Base en1b = r[en - 1];
  const Base enb = r[en];
  energy += m->HairpinInitiation(length) + m->AuGuPenalty(stb, enb);
  const bool all_c = hairpin[st + 1].num_c >= length;

  if (length == 3) {
    if (all_c) energy += m->hairpin_c3_loop;
    return energy;
  }
  energy += m->terminal[r[st]][st1b][en1b][r[en]];

  if ((st1b == U && en1b == U) || (st1b == G && en1b == A))
    energy += m->hairpin_uu_ga_first_mismatch;
  if (st1b == G && en1b == G) energy += m->hairpin_gg_first_mismatch;
  if (all_c) energy += m->hairpin_all_c_a * length + m->hairpin_all_c_b;
  if (stb == G && enb == U && st >= 2 && r[st - 1] == G && r[st - 2] == G)
    energy += m->hairpin_special_gu_closure;

  return energy;
}

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_ENERGY_PRECOMP_H_
