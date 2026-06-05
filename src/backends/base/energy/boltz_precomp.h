// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_BASE_ENERGY_BOLTZ_PRECOMP_H_
#define BACKENDS_BASE_ENERGY_BOLTZ_PRECOMP_H_
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <type_traits>
#include <utility>

#include "api/energy/energy_cfg.h"
#include "backends/base/energy/boltz_model.h"
#include "backends/common/base/boltz_precomp_base.h"
#include "model/energy.h"
#include "model/primary.h"

namespace mrna::md::base {

template <erg::EnergyCfg CFG, class PseudofreeCfgT, class BoltzModelT = BoltzModel>
struct BoltzPrecomp : public BoltzPrecompBase<BoltzModelT> {
  using BoltzPseudofreeCfgT = std::decay_t<decltype(Boltz(std::declval<const PseudofreeCfgT&>()))>;

  PseudofreeCfgT pf;
  BoltzPseudofreeCfgT bpf;

  BoltzPrecomp(Primary r, typename BoltzModelT::Ptr bm, PseudofreeCfgT pf);

  [[nodiscard]] BoltzEnergy Hairpin(int st, int en) const;

  [[nodiscard]] BoltzEnergy TwoLoop(int ost, int oen, int ist, int ien) const;
};

template <erg::EnergyCfg CFG, class PseudofreeCfgT, class BoltzModelT>
BoltzPrecomp<CFG, PseudofreeCfgT, BoltzModelT>::BoltzPrecomp(
    Primary r, typename BoltzModelT::Ptr bm, PseudofreeCfgT pf_)
    : BoltzPrecompBase<BoltzModelT>(std::move(r), std::move(bm)), pf(std::move(pf_)),
      bpf(Boltz(pf)) {
  bpf.Verify(this->r_);
}

template <erg::EnergyCfg CFG, class PseudofreeCfgT, class BoltzModelT>
BoltzEnergy BoltzPrecomp<CFG, PseudofreeCfgT, BoltzModelT>::Hairpin(int st, int en) const {
  const auto& r = this->r_;
  const auto& bm = this->bm_;
  const auto& hairpin = this->hairpin;
  const auto& m = bm->m();
  const int length = en - st - 1;
  assert(length >= HAIRPIN_MIN_SZ);

  BoltzEnergy energy = bpf.UnpairedProd(st + 1, en - 1) * bpf.Paired(st, en);

  if (length <= MAX_SPECIAL_HAIRPIN_SZ && hairpin[st].special[length] > -1)
    return energy * hairpin[st].special[length];

  const Base stb = r[st];
  const Base st1b = r[st + 1];
  const Base en1b = r[en - 1];
  const Base enb = r[en];
  energy *= Boltz(m.HairpinInitiation(length) + m.AuGuPenalty(stb, enb));
  const bool all_c = hairpin[st + 1].num_c >= length;

  if (length == 3) {
    if (all_c) energy *= bm->hairpin_c3_loop;
    return energy;
  }
  energy *= bm->terminal[r[st]][st1b][en1b][r[en]];

  if ((st1b == U && en1b == U) || (st1b == G && en1b == A))
    energy *= bm->hairpin_uu_ga_first_mismatch;
  if (st1b == G && en1b == G) energy *= bm->hairpin_gg_first_mismatch;
  if (all_c) energy *= Boltz(m.hairpin_all_c_a * length) * bm->hairpin_all_c_b;
  if (stb == G && enb == U && st >= 2 && r[st - 1] == G && r[st - 2] == G)
    energy *= bm->hairpin_special_gu_closure;

  return energy;
}

template <erg::EnergyCfg CFG, class PseudofreeCfgT, class BoltzModelT>
BoltzEnergy BoltzPrecomp<CFG, PseudofreeCfgT, BoltzModelT>::TwoLoop(
    int ost, int oen, int ist, int ien) const {
  const auto& r = this->r_;
  const auto& bm = this->bm_;
  const auto& m = bm->m();
  const int toplen = ist - ost - 1;
  const int botlen = oen - ien - 1;
  if (toplen == 0 && botlen == 0)
    return bm->stack[r[ost]][r[ist]][r[ien]][r[oen]] * bpf.Paired(ost, oen);
  if (toplen == 0 || botlen == 0) return bm->template Bulge<CFG>(r, pf, ost, oen, ist, ien);

  BoltzEnergy energy = bm->AuGuPenalty(r[ost], r[oen]) * bm->AuGuPenalty(r[ist], r[ien]) *
      bpf.Paired(ost, oen) * bpf.UnpairedProd(ost + 1, ist - 1) *
      bpf.UnpairedProd(ien + 1, oen - 1);

  if (toplen == 1 && botlen == 1)
    return energy * bm->internal_1x1[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[oen]];
  if (toplen == 1 && botlen == 2)
    return energy *
        bm->internal_1x2[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]][r[oen]];
  if (toplen == 2 && botlen == 1)
    return energy *
        bm->internal_1x2[r[ien]][r[ien + 1]][r[oen]][r[ost]][r[ost + 1]][r[ost + 2]][r[ist]];
  if (toplen == 2 && botlen == 2)
    return energy *
        bm->internal_2x2[r[ost]][r[ost + 1]][r[ost + 2]][r[ist]][r[ien]][r[ien + 1]]
                        [r[ien + 2]][r[oen]];

  static_assert(TWOLOOP_MAX_SZ <= std::remove_reference_t<decltype(m)>::INITIATION_CACHE_SZ,
      "initiation cache not large enough");
  assert(toplen + botlen < std::remove_reference_t<decltype(m)>::INITIATION_CACHE_SZ);
  energy *= bm->internal_init[toplen + botlen] *
      Boltz(std::min(std::abs(toplen - botlen) * m.internal_asym, NINIO_MAX_ASYM));

  energy *= bm->InternalLoopAuGuPenalty(r[ost], r[oen]);
  energy *= bm->InternalLoopAuGuPenalty(r[ist], r[ien]);

  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2))
    energy *= bm->internal_2x3_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] *
        bm->internal_2x3_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];
  else if (toplen != 1 && botlen != 1)
    energy *= bm->internal_other_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] *
        bm->internal_other_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];

  return energy;
}

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_ENERGY_BOLTZ_PRECOMP_H_
