// Copyright 2016 Eliot Courtney.
#include "backends/base/energy/model.h"

#include <cassert>
#include <deque>
#include <memory>
#include <string>
#include <utility>

#include "api/ctx/backend_cfg.h"
#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "backends/common/base/model_base.h"
#include "backends/common/static.h"
#include "model/base.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "model/structure.h"

namespace mrna::md::base {

bool Model::IsSupported(const BackendCfg& backend_cfg, const erg::EnergyCfg& /*cfg*/,
    const erg::PseudofreeCfg& /*pf*/, std::string* reason) {
  if (backend_cfg.energy_model == erg::EnergyModelKind::T22) {
    if (reason) *reason = "T22 energy model not supported by base backend";
    return false;
  }
  return true;
}

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
Energy Model::SubEnergyInternal(const Primary& r, const Secondary& s, const PseudofreeCfgT& pf,
    int st, int en, bool use_given_ctds, Ctds* ctd, std::unique_ptr<Structure>* struc) const {
  assert(en >= st);
  const bool exterior_loop = s[st] != en;
  Energy energy = ZERO_E;

  // Look for branches inside.
  std::deque<int> branches;
  for (int i = st; i <= en; ++i) {
    const int pair = s[i];
    assert(s.IsUnpaired(i) || (pair <= en && s[pair] == i));
    if ((i != st || pair != en) && (i != en || pair != st) && s.IsPaired(i)) {
      branches.push_back(i);
      // Skip ahead.
      i = pair;
    }
  }

  if (exterior_loop || branches.size() >= 2) {
    // Multiloop.
    energy += MultiloopEnergy<CFG>(r, s, pf, st, en, &branches, use_given_ctds, ctd, struc);
  } else if (branches.empty()) {
    // Hairpin loop.
    assert(en - st - 1 >= 3);
    energy += Hairpin<CFG>(r, pf, st, en, struc);
  } else {
    assert(branches.size() == 1);
    const int loop_st = branches.front();
    const int loop_en = s[branches.front()];
    energy += TwoLoop<CFG>(r, pf, st, en, loop_st, loop_en, struc);
  }

  if (struc) (*struc)->set_self_energy(energy);
  // Add energy from children.
  for (auto i : branches) {
    if (struc) {
      std::unique_ptr<Structure> sstruc;
      energy += SubEnergyInternal<CFG>(r, s, pf, i, s[i], use_given_ctds, ctd, &sstruc);
      (*struc)->AddBranch(std::move(sstruc));
    } else {
      energy += SubEnergyInternal<CFG>(r, s, pf, i, s[i], use_given_ctds, ctd, nullptr);
    }
  }
  if (struc) (*struc)->set_total_energy(energy);

  return energy;
}

// If (st, en) is not paired, treated as an exterior loop.
// If `ctd` is non-null, use the given ctds.
template <erg::EnergyCfg CFG, class PseudofreeCfgT>
EnergyResult Model::SubEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd,
    const PseudofreeCfgT& pf, int st, int en, bool build_structure) const {
  ModelBase::Verify(r, s, given_ctd);
  pf.Verify(r);

  const bool use_given_ctds = given_ctd;
  auto ctd = use_given_ctds ? Ctds(*given_ctd) : Ctds(r.size());

  std::unique_ptr<Structure> struc;
  auto energy = SubEnergyInternal<CFG>(
      r, s, pf, st, en, use_given_ctds, &ctd, build_structure ? &struc : nullptr);
  return {energy, std::move(ctd), std::move(struc)};
}

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
EnergyResult Model::TotalEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd,
    const PseudofreeCfgT& pf, bool build_structure) const {
  const int N = r.size();
  auto res = SubEnergy<CFG>(r, s, given_ctd, pf, 0, N - 1, build_structure);
  if (!given_ctd && s[0] == N - 1) res.ctd[0] = CTD_UNUSED;
  if (s[0] == N - 1 && IsAuPair(r[0], r[s[0]])) {
    res.energy += au_penalty;
    if (res.struc) {
      res.struc->AddNote("{}e - top level AU penalty", au_penalty);
      res.struc->set_self_energy(res.struc->self_energy() + au_penalty);
      res.struc->set_total_energy(res.struc->total_energy() + au_penalty);
    }
  }
  if (s[0] == N - 1 && IsGuPair(r[0], r[s[0]])) {
    res.energy += gu_penalty;
    if (res.struc) {
      res.struc->AddNote("{}e - top level GU penalty", gu_penalty);
      res.struc->set_self_energy(res.struc->self_energy() + gu_penalty);
      res.struc->set_total_energy(res.struc->total_energy() + gu_penalty);
    }
  }
  return res;
}

#define MRNA_INSTANTIATE(CFG, PF)                                                      \
  template EnergyResult Model::SubEnergy<CFG, PF>(                                     \
      const Primary&, const Secondary&, const Ctds*, const PF&, int, int, bool) const; \
  template EnergyResult Model::TotalEnergy<CFG, PF>(                                   \
      const Primary&, const Secondary&, const Ctds*, const PF&, bool) const;
MRNA_STATIC_FOR_EACH_ENERGY_CFG_PF(MRNA_INSTANTIATE)
#undef MRNA_INSTANTIATE

}  // namespace mrna::md::base
