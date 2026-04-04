// Copyright 2021 Eliot Courtney.
#include "backends/brute/brute.h"

#include <fmt/core.h>

#include <iterator>
#include <string>
#include <utility>

#include "api/ctx/algorithm.h"
#include "api/energy/energy.h"
#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "api/trace/trace.h"
#include "model/branch.h"
#include "model/constants.h"
#include "model/energy.h"
#include "util/error.h"
#include "util/log.h"

namespace mrna::md::brute {

Brute::Brute(const Primary& r, BackendModelPtr m, BackendCfg backend_cfg, erg::EnergyCfg cfg,
    erg::PseudofreeCfg pf, BruteCfg brute_cfg)
    : r_(r), m_(std::move(m)), bm_(Boltz(m_)), underlying_(Underlying(bm_)),
      backend_cfg_(std::move(backend_cfg)), pf_(std::move(pf)), energy_cfg_(cfg),
      pfn_energy_cfg_(cfg), brute_cfg_(brute_cfg), s_(r_.size()), ctd_(r_.size()) {
  std::string reason;
  verify(BackendIsSupported(GetBackendKind(m_), backend_cfg_, energy_cfg_, pf_, &reason),
      "backend {} does not support configuration {}: {}", GetBackendKind(m_), energy_cfg_, reason);
  pfn_energy_cfg_.bulge_states = false;

  verify(brute_cfg_.subopt_cfg.time_secs < 0,
      "brute force does not support time limit for suboptimal folding");
}

BruteResult Brute::Run() {
  // Preconditions:
  static_assert(CTD_SIZE < (1 << CTD_MAX_BITS), "need increase ctd bits for brute force");

  logdebug("brute {} with {}, {}, {}", funcname(), energy_cfg_, brute_cfg_, pf_);

  if (brute_cfg_.pfn) {
    // Plus one to N, since the packed representation reserves an all-ones sentinel.
    verify(r_.size() + 1 < (1 << PT_MAX_BITS), "sequence too long for brute force partition");
    res_.pfn.q = 0;
    res_.pfn.p = BoltzSums(r_.size(), 0);
  }
  // Add base pairs in order of increasing st, then en.
  for (int st = 0; st < r_.size(); ++st) {
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < r_.size(); ++en) {
      if (CanPair(energy_cfg_, m_, r_, st, en)) pairs_.emplace_back(st, en);
    }
  }
  Dfs(0);

  if (brute_cfg_.pfn) res_.pfn.RecomputeProb();

  return std::move(res_);
}

void Brute::Dfs(int idx) {
  if (idx == int(pairs_.size())) {
    // Precompute whether things are multiloops or not.
    branch_count_ = GetBranchCounts(s_);
    AddAllCombinations(0);
    return;
  }
  // Don't take this base pair.
  Dfs(idx + 1);

  // Take this base pair.
  bool can_take = true;
  const auto& p = pairs_[idx];
  // Only need to check in the range of this base pair. Since we ordered by
  // increasing st, anything at or after this will either be the start of something starting at st,
  // or something ending, both of which conflict with this base pair.
  for (int i = p.st; i <= p.en; ++i) {
    if (s_.IsPaired(i)) {
      can_take = false;
      break;
    }
  }
  if (can_take) {
    s_[p.st] = p.en;
    s_[p.en] = p.st;
    Dfs(idx + 1);
    s_[p.st] = INVALID_INDEX;
    s_[p.en] = INVALID_INDEX;
  }
}

void Brute::AddAllCombinations(int idx) {
  const int N = r_.size();
  // Base case
  if (idx == N) {
    if (brute_cfg_.pfn) {
      auto energy = TotalEnergy(underlying_, r_, s_, &ctd_, pfn_energy_cfg_, pf_).energy;
      res_.pfn.q += energy.Boltz();
      for (int i = 0; i < N; ++i) {
        if (s_.IsOpeningPair(i)) {
          const auto inside_structure = BuildInsideStructure(i, s_[i], N);
          const auto outside_structure = BuildOutsideStructure(i, s_[i], N);
          const bool inside_new = !substructure_map_.Find(inside_structure);
          const bool outside_new = !substructure_map_.Find(outside_structure);
          if (inside_new || outside_new) {
            const Energy inside_energy =
                SubEnergy(underlying_, r_, s_, &ctd_, pfn_energy_cfg_, pf_, i, s_[i]).energy;
            if (inside_new) {
              res_.pfn.p[i][s_[i]] += inside_energy.Boltz();
              substructure_map_.Insert(inside_structure, Nothing());
            }
            if (outside_new) {
              res_.pfn.p[s_[i]][i] += (energy - inside_energy).Boltz();
              substructure_map_.Insert(outside_structure, Nothing());
            }
          }
        }
      }
    }
    if (brute_cfg_.subopt) {
      auto energy = TotalEnergy(m_, r_, s_, &ctd_, energy_cfg_, pf_).energy;
      PruneInsertSubopt(energy);
    }
    return;
  }

  // If we already set this, this isn't a valid base pair, it's not part of a multiloop, can't set
  // ctds so continue.
  const int st = idx;
  if (ctd_[idx] != CTD_NA || s_.IsUnpaired(st) || branch_count_[idx] < 2) {
    AddAllCombinations(idx + 1);
    return;
  }

  const int en = s_[st];
  const bool lu_exists = st != 0 && s_.IsUnpaired(st - 1);
  const bool lu_shared = lu_exists && st > 1 && s_.IsPaired(st - 2);
  const bool lu_usable = lu_exists &&
      (!lu_shared ||
          (ctd_[s_[st - 2]] != CTD_3_DANGLE && ctd_[s_[st - 2]] != CTD_MISMATCH &&
              ctd_[s_[st - 2]] != CTD_RC_WITH_PREV));
  const bool ru_exists = en + 1 < N && s_.IsUnpaired(en + 1);
  const bool ru_shared = ru_exists && en + 2 < N && s_.IsPaired(en + 2);
  const bool ru_usable = ru_exists &&
      (!ru_shared ||
          (ctd_[en + 2] != CTD_5_DANGLE && ctd_[en + 2] != CTD_MISMATCH &&
              ctd_[en + 2] != CTD_LCOAX_WITH_NEXT));
  // Even if the next branch is an outer branch, everything will be magically handled.

  // CTD_NONE or D2 CTD;
  Ctd none_ctd = CTD_UNUSED;
  if (energy_cfg_.UseD2()) {
    const bool lspace = st != 0;
    const bool rspace = en + 1 < N;
    if (lspace && rspace) {
      none_ctd = CTD_MISMATCH;
    } else if (rspace) {
      none_ctd = CTD_3_DANGLE;
    } else if (lspace) {
      none_ctd = CTD_5_DANGLE;
    }
  }
  ctd_[idx] = none_ctd;
  AddAllCombinations(idx + 1);

  if (energy_cfg_.UseDangleMismatch()) {
    // CTD_3_DANGLE
    if (ru_usable) {
      ctd_[idx] = CTD_3_DANGLE;
      AddAllCombinations(idx + 1);
    }

    // CTD_5_DANGLE
    if (lu_usable) {
      ctd_[idx] = CTD_5_DANGLE;
      AddAllCombinations(idx + 1);
    }

    // CTD_MISMATCH
    if (ru_usable && lu_usable) {
      ctd_[idx] = CTD_MISMATCH;
      AddAllCombinations(idx + 1);
    }
  }

  if (energy_cfg_.UseCoaxialStacking()) {
    // Check that the next branch hasn't been set already. If it's unused or na, try re-writing it.
    // CTD_LCOAX_WITH_NEXT
    if (lu_usable && ru_usable && ru_shared) {
      auto prevval = ctd_[en + 2];
      if (prevval == CTD_UNUSED || prevval == CTD_NA) {
        ctd_[idx] = CTD_LCOAX_WITH_NEXT;
        ctd_[en + 2] = CTD_LCOAX_WITH_PREV;
        AddAllCombinations(idx + 1);
        ctd_[en + 2] = prevval;
      }
    }

    // Check that the previous branch hasn't been set already.
    // CTD_RC_WITH_PREV
    if (lu_usable && lu_shared && ru_usable) {
      auto prevval = ctd_[s_[idx - 2]];
      if (prevval == CTD_UNUSED || prevval == CTD_NA) {
        ctd_[idx] = CTD_RC_WITH_PREV;
        ctd_[s_[idx - 2]] = CTD_RC_WITH_NEXT;
        AddAllCombinations(idx + 1);
        ctd_[s_[idx - 2]] = prevval;
      }
    }

    // CTD_FCOAX_WITH_NEXT
    if (en + 1 < N && s_.IsPaired(en + 1)) {
      auto prevval = ctd_[en + 1];
      if (prevval == CTD_UNUSED || prevval == CTD_NA) {
        ctd_[idx] = CTD_FCOAX_WITH_NEXT;
        ctd_[en + 1] = CTD_FCOAX_WITH_PREV;
        AddAllCombinations(idx + 1);
        ctd_[en + 1] = prevval;
      }
    }
  }

  // Reset back to NA.
  ctd_[idx] = CTD_NA;
}

void Brute::PruneInsertSubopt(Energy e) {
  const bool has_room = static_cast<int64_t>(res_.subopts.size()) < brute_cfg_.subopt_cfg.strucs;
  const bool is_better = res_.subopts.empty() || res_.subopts.rbegin()->energy > e;
  if (has_room || is_better)
    res_.subopts.insert(subopt::SuboptResult(e, trace::TraceResult(Secondary(s_), Ctds(ctd_))));

  // Prune values that exceed the number of structures:
  if (static_cast<int64_t>(res_.subopts.size()) > brute_cfg_.subopt_cfg.strucs)
    res_.subopts.erase(--res_.subopts.end());

  // Prune values that exceed the delta:
  while (!res_.subopts.empty()) {
    if (res_.subopts.rbegin()->energy - res_.subopts.begin()->energy <= brute_cfg_.subopt_cfg.delta)
      break;
    res_.subopts.erase(--res_.subopts.end());
  }
}

Brute::SubstructureId Brute::WriteBits(int st, int en, int N, bool inside) {
  static_assert(PT_MAX_BITS + CTD_MAX_BITS <= 16, "substructure block does not fit in uint16_t");
  // Zero initialise. Unpaired is all-ones in the packed representation once PT_MASK is applied.
  SubstructureId struc = {};
  uint32_t b = 0;
  for (int i = 0; i < N; ++i, b += PT_MAX_BITS + CTD_MAX_BITS) {
    if (inside && (i < st || i > en)) continue;
    if (!inside && i > st && i < en) continue;

    // Don't include the ctd value at st, since that's for the outside.
    // Don't include the ctd value at en, since that's for the inside.
    auto ctd = ctd_[i];
    if ((inside && i == st) || (!inside && i == en)) ctd = CTD_NA;
    auto pack = uint16_t((uint32_t(s_[i]) & PT_MASK) << CTD_MAX_BITS | (ctd & CTD_MASK));

    const uint32_t chunk = b / 16;
    const uint32_t bit = b & 15;
    struc.bits[chunk] |= uint16_t(pack << bit);
    const uint32_t space = 16 - bit;
    if (space < CTD_MAX_BITS + PT_MAX_BITS) struc.bits[chunk + 1] = uint16_t(pack >> space);
  }
  return struc;
}

Brute::SubstructureId Brute::BuildInsideStructure(int st, int en, int N) {
  return WriteBits(st, en, N, /*inside=*/true);
}

Brute::SubstructureId Brute::BuildOutsideStructure(int st, int en, int N) {
  return WriteBits(st, en, N, /*inside=*/false);
}

}  // namespace mrna::md::brute
