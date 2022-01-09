// Copyright 2021 Eliot Courtney.
#include "compute/brute.h"

#include "compute/constants.h"
#include "util/macros.h"

namespace mrna {

namespace internal {

std::vector<int> GetBranchCounts(const std::vector<int>& p) {
  std::vector<int> branch_count(p.size(), 0);
  std::stack<int> q;
  for (int i = 0; i < static_cast<int>(p.size()); ++i) {
    if (p[i] == -1) continue;
    if (p[i] > i) {
      // Exterior loop counts a multiloop for CTDs.
      if (q.empty()) branch_count[i] = 2;
      q.push(i);

      // Look at all the children.
      int count = 0;
      for (int j = i + 1; j < p[i]; ++j) {
        if (p[j] != -1) {
          j = p[j];
          count++;
        }
      }
      for (int j = i + 1; j < p[i]; ++j) {
        if (p[j] != -1) {
          branch_count[j] = count;
          j = p[j];
        }
      }
      branch_count[p[i]] = count;
    } else {
      q.pop();
    }
  }
  return branch_count;
}

}  // namespace internal

BruteForce::SubstructureId BruteForce::WriteBits(int st, int en, int N, bool inside) {
  static_assert(PT_MAX_BITS + CTD_MAX_BITS <= 16, "substructure block does not fit in uint16_t");
  static_assert((-1 & PT_MASK) == PT_MASK, "mfw not a two's complement machine");
  SubstructureId struc = {};  // Zero initialise.
  for (int i = 0, b = 0; i < N; ++i, b += PT_MAX_BITS + CTD_MAX_BITS) {
    if (inside && (i < st || i > en)) continue;
    if (!inside && i > st && i < en) continue;
    uint16_t pack = uint16_t((p_[i] & PT_MASK) << CTD_MAX_BITS | (ctd_[i] & CTD_MASK));

    int byte = b / 16;
    int bit = b & 15;
    struc.bits[byte] = uint16_t(struc.bits[byte] | (pack << bit));
    int space = 16 - bit;
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

void BruteForce::AddAllCombinations(int idx) {
  const int N = static_cast<int>(r_.size());
  // Base case
  if (idx == N) {
    auto computed = energy::ComputeEnergyWithCtds({{r_, p_}, ctd_, 0}, em_);
    if (res_.compute_partition) {
      BoltzEnergy boltz = Boltz(computed.energy);
      res_.partition.q += boltz;
      for (int i = 0; i < N; ++i) {
        if (i < p_[i]) {
          const auto inside_structure = BuildInsideStructure(i, p_[i], N);
          const auto outside_structure = BuildOutsideStructure(i, p_[i], N);
          const bool inside_new = !substructure_map_.Find(inside_structure);
          const bool outside_new = !substructure_map_.Find(outside_structure);
          if (inside_new || outside_new) {
            Energy inside_energy = energy::ComputeSubstructureEnergy(
                computed, false, i, p_[i], em_);  // TODO optimisation?
            if (inside_new) {
              res_.partition.p[i][p_[i]][0] += Boltz(inside_energy);
              substructure_map_.Insert(inside_structure, Nothing());
            }
            if (outside_new) {
              res_.partition.p[p_[i]][i][0] += Boltz(computed.energy - inside_energy);
              substructure_map_.Insert(outside_structure, Nothing());
            }
          }
          res_.probabilities[i][p_[i]][0] += boltz;
        }
      }
    } else {
      if (static_cast<int>(res_.best_computeds.size()) < res_.max_structures ||
          res_.best_computeds.rbegin()->energy > computed.energy)
        res_.best_computeds.insert(std::move(computed));
      if (static_cast<int>(res_.best_computeds.size()) > res_.max_structures)
        res_.best_computeds.erase(--res_.best_computeds.end());
    }
    return;
  }

  // If we already set this, this isn't a valid base pair, it's not part of a multiloop, can't set
  // ctds so continue.
  if (ctd_[idx] != CTD_NA || p_[idx] == -1 || res_.branch_count[idx] < 2) {
    AddAllCombinations(idx + 1);
    return;
  }

  const bool lu_exists = idx - 1 >= 0 && p_[idx - 1] == -1;
  const bool lu_shared = lu_exists && idx - 2 >= 0 && p_[idx - 2] != -1;
  const bool lu_usable = lu_exists &&
      (!lu_shared ||
          (ctd_[p_[idx - 2]] != CTD_3_DANGLE && ctd_[p_[idx - 2]] != CTD_MISMATCH &&
              ctd_[p_[idx - 2]] != CTD_RCOAX_WITH_PREV));
  const bool ru_exists = p_[idx] + 1 < N && p_[p_[idx] + 1] == -1;
  const bool ru_shared = ru_exists && p_[idx] + 2 < N && p_[p_[idx] + 2] != -1;
  const bool ru_usable = ru_exists &&
      (!ru_shared ||
          (ctd_[p_[idx] + 2] != CTD_5_DANGLE && ctd_[p_[idx] + 2] != CTD_MISMATCH &&
              ctd_[p_[idx] + 2] != CTD_LCOAX_WITH_NEXT));
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
    auto prevval = ctd_[p_[idx] + 2];
    if (prevval == CTD_UNUSED || prevval == CTD_NA) {
      ctd_[idx] = CTD_LCOAX_WITH_NEXT;
      ctd_[p_[idx] + 2] = CTD_LCOAX_WITH_PREV;
      AddAllCombinations(idx + 1);
      ctd_[p_[idx] + 2] = prevval;
    }
  }

  // Check that the previous branch hasn't been set already.
  // CTD_RCOAX_WITH_PREV
  if (lu_usable && lu_shared && ru_usable) {
    auto prevval = ctd_[p_[idx - 2]];
    if (prevval == CTD_UNUSED || prevval == CTD_NA) {
      ctd_[idx] = CTD_RCOAX_WITH_PREV;
      ctd_[p_[idx - 2]] = CTD_RCOAX_WITH_NEXT;
      AddAllCombinations(idx + 1);
      ctd_[p_[idx - 2]] = prevval;
    }
  }

  // CTD_FCOAX_WITH_NEXT
  if (p_[idx] + 1 < N && p_[p_[idx] + 1] != -1) {
    auto prevval = ctd_[p_[idx] + 1];
    if (prevval == CTD_UNUSED || prevval == CTD_NA) {
      ctd_[idx] = CTD_FCOAX_WITH_NEXT;
      ctd_[p_[idx] + 1] = CTD_FCOAX_WITH_PREV;
      AddAllCombinations(idx + 1);
      ctd_[p_[idx] + 1] = prevval;
    }
  }

  // Reset back to NA.
  ctd_[idx] = CTD_NA;
}

void BruteForce::Dfs(int idx) {
  if (idx == static_cast<int>(res_.base_pairs.size())) {
    // Small optimisation for case when we're just getting one structure.
    if (res_.max_structures == 1 && !res_.compute_partition) {
      auto computed = energy::ComputeEnergy({r_, p_}, em_);
      if (res_.best_computeds.empty() || computed.energy < res_.best_computeds.begin()->energy)
        res_.best_computeds.insert(std::move(computed));
      if (res_.best_computeds.size() == 2) res_.best_computeds.erase(--res_.best_computeds.end());
    } else {
      // Precompute whether things are multiloops or not.
      res_.branch_count = internal::GetBranchCounts(p_);
      AddAllCombinations(0);
    }

    return;
  }
  // Don't take this base pair.
  Dfs(idx + 1);

  // Take this base pair.
  bool can_take = true;
  const auto& pair = res_.base_pairs[idx];
  // Only need to check in the range of this base pair. Since we ordered by
  // increasing st, anything at or after this will either be the start of something starting at st,
  // or something ending, both of which conflict with this base pair.
  for (int i = pair.first; i <= pair.second; ++i) {
    if (p_[i] != -1) {
      can_take = false;
      break;
    }
  }
  if (can_take) {
    p_[pair.first] = pair.second;
    p_[pair.second] = pair.first;
    Dfs(idx + 1);
    p_[pair.first] = -1;
    p_[pair.second] = -1;
  }
}

BruteForce::Result BruteForce::Run(const Primary& r, const energy::EnergyModel& em,
    int max_structures, bool compute_partition, bool allow_lonely_pairs) {
  // Preconditions:
  static_assert(CTD_SIZE < (1 << CTD_MAX_BITS), "need increase ctd bits for brute force");

  r_ = r;
  const int N = static_cast<int>(r.size());

  em_ = em;
  // TODO: improve this - see similar code in  subopt1
  p_.resize(N);
  std::fill(p_.begin(), p_.end(), -1);
  ctd_.resize(N);
  std::fill(ctd_.begin(), ctd_.end(), CTD_NA);

  // TODO: Cleanup here
  res_.max_structures = max_structures == -1 ? MAX_STRUCTURES : max_structures;
  res_.compute_partition = compute_partition;
  if (res_.compute_partition) {
    // Plus one to N, since -1 takes up a spot.
    verify(N + 1 < (1 << PT_MAX_BITS), "sequence too long for brute force partition");
    res_.partition.q = 0;
    res_.partition.p = Array3D<BoltzEnergy, 1>(r.size(), 0);
    res_.probabilities = Array3D<BoltzEnergy, 1>(r.size(), 0);
  }
  // Add base pairs in order of increasing st, then en.
  for (int st = 0; st < N; ++st) {
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < N; ++en) {
      bool allowed = allow_lonely_pairs ? CanPair(r[st], r[en]) : ViableFoldingPair(r, st, en);
      if (allowed) res_.base_pairs.emplace_back(st, en);
    }
  }
  Dfs(0);

  if (res_.compute_partition) {
    // Fill probabilities from partition function.
    for (int i = 0; i < N; ++i)
      for (int j = i; j < N; ++j) res_.probabilities[i][j][0] /= res_.partition.q;
  }

  return std::move(res_);
}

}  // namespace mrna
