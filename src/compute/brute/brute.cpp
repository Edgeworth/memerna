// Copyright 2021 Eliot Courtney.
#include "compute/brute/brute.h"

#include <iterator>
#include <memory>
#include <utility>

#include "compute/energy/branch.h"
#include "compute/energy/energy.h"
#include "compute/partition/partition.h"
#include "compute/subopt/config.h"
#include "compute/subopt/subopt.h"
#include "compute/traceback/traceback.h"
#include "model/model.h"
#include "util/error.h"

namespace mrna::brute {

BruteForce::BruteForce(const Primary& r, energy::EnergyModelPtr em, BruteCfg cfg)
    : r_(r), em_(std::move(em)), cfg_(cfg), s_(r_.size()), ctd_(r_.size()) {}

BruteResult BruteForce::Run() {
  // Preconditions:
  static_assert(CTD_SIZE < (1 << CTD_MAX_BITS), "need increase ctd bits for brute force");

  if (cfg_.part) {
    // Plus one to N, since -1 takes up a spot.
    verify(r_.size() + 1 < (1 << PT_MAX_BITS), "sequence too long for brute force partition");
    res_.part.q = 0;
    res_.part.p = BoltzSums(r_.size(), 0);
  }
  // Add base pairs in order of increasing st, then en.
  for (int st = 0; st < static_cast<int>(r_.size()); ++st) {
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < static_cast<int>(r_.size()); ++en) {
      if (em_->CanPair(r_, st, en)) pairs_.emplace_back(st, en);
    }
  }
  Dfs(0);

  if (cfg_.part) res_.prob = res_.part.Prob();

  return std::move(res_);
}

void BruteForce::Dfs(int idx) {
  if (idx == static_cast<int>(pairs_.size())) {
    // Precompute whether things are multiloops or not.
    branch_count_ = energy::GetBranchCounts(s_);
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
  for (int i = p.first; i <= p.second; ++i) {
    if (s_[i] != -1) {
      can_take = false;
      break;
    }
  }
  if (can_take) {
    s_[p.first] = p.second;
    s_[p.second] = p.first;
    Dfs(idx + 1);
    s_[p.first] = -1;
    s_[p.second] = -1;
  }
}

void BruteForce::AddAllCombinations(int idx) {
  const int N = static_cast<int>(r_.size());
  // Base case
  if (idx == N) {
    auto energy = em_->TotalEnergy(r_, s_, &ctd_).energy;
    if (cfg_.part) {
      BoltzEnergy boltz = Boltz(energy);
      res_.part.q += boltz;
      for (int i = 0; i < N; ++i) {
        if (i < s_[i]) {
          const auto inside_structure = BuildInsideStructure(i, s_[i], N);
          const auto outside_structure = BuildOutsideStructure(i, s_[i], N);
          const bool inside_new = !substructure_map_.Find(inside_structure);
          const bool outside_new = !substructure_map_.Find(outside_structure);
          if (inside_new || outside_new) {
            Energy inside_energy = em_->SubstructureEnergy(r_, s_, &ctd_, i, s_[i]).energy;
            if (inside_new) {
              res_.part.p[i][s_[i]] += Boltz(inside_energy);
              substructure_map_.Insert(inside_structure, Nothing());
            }
            if (outside_new) {
              res_.part.p[s_[i]][i] += Boltz(energy - inside_energy);
              substructure_map_.Insert(outside_structure, Nothing());
            }
          }
        }
      }
    }
    if (cfg_.subopt) PruneInsertSubopt(energy);
    return;
  }

  // If we already set this, this isn't a valid base pair, it's not part of a multiloop, can't set
  // ctds so continue.
  if (ctd_[idx] != CTD_NA || s_[idx] == -1 || branch_count_[idx] < 2) {
    AddAllCombinations(idx + 1);
    return;
  }

  const bool lu_exists = idx - 1 >= 0 && s_[idx - 1] == -1;
  const bool lu_shared = lu_exists && idx - 2 >= 0 && s_[idx - 2] != -1;
  const bool lu_usable = lu_exists &&
      (!lu_shared ||
          (ctd_[s_[idx - 2]] != CTD_3_DANGLE && ctd_[s_[idx - 2]] != CTD_MISMATCH &&
              ctd_[s_[idx - 2]] != CTD_RCOAX_WITH_PREV));
  const bool ru_exists = s_[idx] + 1 < N && s_[s_[idx] + 1] == -1;
  const bool ru_shared = ru_exists && s_[idx] + 2 < N && s_[s_[idx] + 2] != -1;
  const bool ru_usable = ru_exists &&
      (!ru_shared ||
          (ctd_[s_[idx] + 2] != CTD_5_DANGLE && ctd_[s_[idx] + 2] != CTD_MISMATCH &&
              ctd_[s_[idx] + 2] != CTD_LCOAX_WITH_NEXT));
  // Even if the next branch is an outer branch, everything will be magically handled.
  // CTD_UNUSED
  ctd_[idx] = CTD_UNUSED;
  AddAllCombinations(idx + 1);

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

  // Check that the next branch hasn't been set already. If it's unused or na, try re-writing it.
  // CTD_LCOAX_WITH_NEXT
  if (lu_usable && ru_usable && ru_shared) {
    auto prevval = ctd_[s_[idx] + 2];
    if (prevval == CTD_UNUSED || prevval == CTD_NA) {
      ctd_[idx] = CTD_LCOAX_WITH_NEXT;
      ctd_[s_[idx] + 2] = CTD_LCOAX_WITH_PREV;
      AddAllCombinations(idx + 1);
      ctd_[s_[idx] + 2] = prevval;
    }
  }

  // Check that the previous branch hasn't been set already.
  // CTD_RCOAX_WITH_PREV
  if (lu_usable && lu_shared && ru_usable) {
    auto prevval = ctd_[s_[idx - 2]];
    if (prevval == CTD_UNUSED || prevval == CTD_NA) {
      ctd_[idx] = CTD_RCOAX_WITH_PREV;
      ctd_[s_[idx - 2]] = CTD_RCOAX_WITH_NEXT;
      AddAllCombinations(idx + 1);
      ctd_[s_[idx - 2]] = prevval;
    }
  }

  // CTD_FCOAX_WITH_NEXT
  if (s_[idx] + 1 < N && s_[s_[idx] + 1] != -1) {
    auto prevval = ctd_[s_[idx] + 1];
    if (prevval == CTD_UNUSED || prevval == CTD_NA) {
      ctd_[idx] = CTD_FCOAX_WITH_NEXT;
      ctd_[s_[idx] + 1] = CTD_FCOAX_WITH_PREV;
      AddAllCombinations(idx + 1);
      ctd_[s_[idx] + 1] = prevval;
    }
  }

  // Reset back to NA.
  ctd_[idx] = CTD_NA;
}

void BruteForce::PruneInsertSubopt(Energy e) {
  bool has_room = static_cast<int>(res_.subopts.size()) < cfg_.subopt_cfg.strucs;
  bool is_better = res_.subopts.empty() || res_.subopts.rbegin()->energy > e;
  if (has_room || is_better)
    res_.subopts.insert(subopt::SuboptResult(e, tb::TracebackResult(Secondary(s_), Ctds(ctd_))));

  // Prune values that exceed the number of structures:
  if (static_cast<int>(res_.subopts.size()) > cfg_.subopt_cfg.strucs)
    res_.subopts.erase(--res_.subopts.end());

  // Prune values that exceed the delta:
  while (!res_.subopts.empty()) {
    if (res_.subopts.rbegin()->energy - res_.subopts.begin()->energy <= cfg_.subopt_cfg.delta)
      break;
    res_.subopts.erase(--res_.subopts.end());
  }
}

BruteForce::SubstructureId BruteForce::WriteBits(int st, int en, int N, bool inside) {
  static_assert(PT_MAX_BITS + CTD_MAX_BITS <= 16, "substructure block does not fit in uint16_t");
  SubstructureId struc = {};  // Zero initialise.
  uint32_t b = 0;
  for (int i = 0; i < N; ++i, b += PT_MAX_BITS + CTD_MAX_BITS) {
    if (inside && (i < st || i > en)) continue;
    if (!inside && i > st && i < en) continue;
    auto pack = uint16_t((uint32_t(s_[i]) & PT_MASK) << CTD_MAX_BITS | (ctd_[i] & CTD_MASK));

    uint32_t byte = b / 16;
    uint32_t bit = b & 15;
    struc.bits[byte] = uint16_t(struc.bits[byte] | (pack << bit));
    uint32_t space = 16 - bit;
    if (space < CTD_MAX_BITS + PT_MAX_BITS)
      struc.bits[byte + 1] = uint16_t(struc.bits[byte + 1] | (pack >> space));
  }
  return struc;
}

BruteForce::SubstructureId BruteForce::BuildInsideStructure(int st, int en, int N) {
  // Don't include the ctd value at st, since that's for the outside.
  return WriteBits(st + 1, en, N, true);
}

BruteForce::SubstructureId BruteForce::BuildOutsideStructure(int st, int en, int N) {
  // Don't include the ctd value at en, since that's for the inside.
  return WriteBits(st, en + 1, N, false);
}

}  // namespace mrna::brute
