// Copyright 2022 Eliot Courtney.
#include "backends/stack/energy/model.h"

#include <fmt/format.h>

#include <cassert>
#include <deque>
#include <memory>
#include <random>
#include <string>
#include <utility>

#include "api/ctx/backend_cfg.h"
#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "backends/common/base/model_base.h"
#include "backends/common/base/parse.h"
#include "backends/common/model_mixin.h"
#include "backends/common/parse.h"
#include "backends/common/static.h"
#include "model/base.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "model/structure.h"
#include "util/file.h"

namespace mrna::md::stack {

bool Model::IsSupported(const BackendCfg& /*backend_cfg*/, const erg::EnergyCfg& cfg,
    const erg::PseudofreeCfg& /*pf*/, std::string* reason) {
  // STACK supports all energy models (T04, T12, T22).
  if (cfg.ctd != erg::EnergyCfg::Ctd::ALL && cfg.ctd != erg::EnergyCfg::Ctd::NO_COAX &&
      cfg.ctd != erg::EnergyCfg::Ctd::NONE) {
    if (reason) *reason = fmt::format("ctd={} not supported by stack backend", cfg.ctd);
    return false;
  }
  return true;
}

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
Energy Model::SubEnergyInternal(const Primary& r, const Secondary& s, const PseudofreeCfgT& pf,
    int st, int en, int stack_st, int stack_en, bool use_given_ctds, Ctds* ctd,
    std::unique_ptr<Structure>* struc) const {
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

  bool continuous = false;
  if (exterior_loop || branches.size() >= 2) {
    // Multiloop.
    energy += MultiloopEnergy<CFG>(r, s, pf, st, en, &branches, use_given_ctds, ctd, struc);

    // Current stack is terminated.
    if (stack_st != -1 && stack_st != st)
      energy += StackPenalty(r, s, stack_st, stack_en, st, en, struc);
    stack_st = -1;
    stack_en = -1;
  } else if (branches.empty()) {
    // Hairpin loop.
    assert(en - st - 1 >= 3);
    energy += Hairpin<CFG>(r, pf, st, en, struc);

    // Current stack is terminated.
    if (stack_st != -1 && stack_st != st)
      energy += StackPenalty(r, s, stack_st, stack_en, st, en, struc);
    stack_st = -1;
    stack_en = -1;
  } else {
    assert(branches.size() == 1);
    const int loop_st = branches.front();
    const int loop_en = s[branches.front()];
    energy += TwoLoop<CFG>(r, pf, st, en, loop_st, loop_en, struc);

    // Current stack is terminated if it's not continuous, otherwise it's extended.
    if (!IsContinuous(st, en, loop_st, loop_en)) {
      if (stack_st != -1 && stack_st != st)
        energy += StackPenalty(r, s, stack_st, stack_en, st, en, struc);
      stack_st = -1;
      stack_en = -1;
    } else {
      continuous = true;
    }
  }

  if (struc) (*struc)->set_self_energy(energy);
  // Add energy from children.
  for (auto i : branches) {
    const int nst = continuous ? stack_st : i;
    const int nen = continuous ? stack_en : s[i];
    if (struc) {
      std::unique_ptr<Structure> sstruc;
      energy += SubEnergyInternal<CFG>(r, s, pf, i, s[i], nst, nen, use_given_ctds, ctd, &sstruc);
      (*struc)->AddBranch(std::move(sstruc));
    } else {
      energy += SubEnergyInternal<CFG>(r, s, pf, i, s[i], nst, nen, use_given_ctds, ctd, nullptr);
    }
  }
  if (struc) (*struc)->set_total_energy(energy);

  return energy;
}

inline Energy Model::StackPenalty(const Primary& r, const Secondary& s, int ost, int oen, int ist,
    int ien, std::unique_ptr<Structure>* struc) const {
  assert(ost >= 0 && oen < r.size());
  assert(ist > 0 && ien >= 0 && ien + 1 < r.size());
  assert(ist > ost && ien < oen);
  // Check for single unpaired bases to treat as continuous.
  const int ost_next = s.IsUnpaired(ost + 1) ? ost + 2 : ost + 1;
  const int oen_prev = s.IsUnpaired(oen - 1) ? oen - 2 : oen - 1;
  assert(ost_next - ost + oen - oen_prev < 4);

  const int ist_prev = s.IsUnpaired(ist - 1) ? ist - 2 : ist - 1;
  const int ien_next = s.IsUnpaired(ien + 1) ? ien + 2 : ien + 1;
  assert(ist - ist_prev + ien_next - ien < 4);

  auto inner = penultimate_stack[r[ist_prev]][r[ist]][r[ien]][r[ien_next]];
  auto outer = penultimate_stack[r[oen_prev]][r[oen]][r[ost]][r[ost_next]];
  if (struc) {
    (*struc)->AddNote("{}e - inner penultimate penalty at ({}, {})", inner, ist, ien);
    (*struc)->AddNote("{}e - outer penultimate penalty at ({}, {})", outer, ost, oen);
  }

  return inner + outer;
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

  const bool exterior_loop = s[st] != en;
  const int stack_st = exterior_loop ? -1 : st;
  const int stack_en = exterior_loop ? -1 : en;

  std::unique_ptr<Structure> struc;
  auto energy = SubEnergyInternal<CFG>(r, s, pf, st, en, stack_st, stack_en, use_given_ctds, &ctd,
      build_structure ? &struc : nullptr);
  return {energy, std::move(ctd), std::move(struc)};
}

template <erg::EnergyCfg CFG, class PseudofreeCfgT>
EnergyResult Model::TotalEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd,
    const PseudofreeCfgT& pf, bool build_structure) const {
  const int N = r.size();
  auto res = SubEnergy<CFG>(r, s, given_ctd, pf, 0, N - 1, build_structure);
  if (!given_ctd && s[0] == N - 1) res.ctd[0] = CTD_UNUSED;
  if (s[0] == N - 1) {
    auto extra_energy = ZERO_E;
    if (IsAuPair(r[0], r[s[0]])) {
      extra_energy += au_penalty;
      if (res.struc) res.struc->AddNote("{}e - top level AU penalty", au_penalty);
    }
    if (IsGuPair(r[0], r[s[0]])) {
      extra_energy += gu_penalty;
      if (res.struc) res.struc->AddNote("{}e - top level GU penalty", gu_penalty);
    }
    res.energy += extra_energy;
    if (res.struc) {
      res.struc->set_self_energy(res.struc->self_energy() + extra_energy);
      res.struc->set_total_energy(res.struc->total_energy() + extra_energy);
    }
  }

  return res;
}

void Model::LoadFromModelPath(const std::string& path) {
  base::LoadFromModelPath(*this, path);
  Parse4Map(OpenFile(path + "/penultimate_stacking.data"), penultimate_stack);
}

void Model::LoadRandom(const BackendCfg& cfg, const RandomModelCfg& random_cfg, std::mt19937& eng) {
  const auto min_energy = random_cfg.min_energy.ToFlt<double>();
  const auto max_energy = random_cfg.max_energy.ToFlt<double>();
  base::LoadRandomModel(
      *this, eng, min_energy, max_energy, RAND_MAX_HAIRPIN_SZ, RAND_MAX_NUM_HAIRPIN);

  // Only randomize penultimate_stack for T22, which uses penultimate stacking.
  if (cfg.energy_model == erg::EnergyModelKind::T22) {
    // penultimate_stack is dependent on the direction, so 180 degree rotations
    // don't have to be the same.
    std::uniform_real_distribution<double> energy_dist(min_energy, max_energy);
    RANDOMISE_DATA((*this), penultimate_stack);
  }
}

#define MRNA_INSTANTIATE(CFG, PF)                                                      \
  template EnergyResult Model::SubEnergy<CFG, PF>(                                     \
      const Primary&, const Secondary&, const Ctds*, const PF&, int, int, bool) const; \
  template EnergyResult Model::TotalEnergy<CFG, PF>(                                   \
      const Primary&, const Secondary&, const Ctds*, const PF&, bool) const;
MRNA_STATIC_FOR_EACH_ENERGY_CFG_PF(MRNA_INSTANTIATE)
#undef MRNA_INSTANTIATE

}  // namespace mrna::md::stack
