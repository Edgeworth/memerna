// Copyright 2016 Eliot Courtney.
#include "fold/brute_fold.h"

#include <set>
#include <stack>
#include <utility>
#include <vector>

#include "energy/energy_globals.h"
#include "energy/structure.h"
#include "fold/fold.h"
#include "fold/fold_globals.h"
#include "parsing.h"
#include "splaymap.h"

namespace memerna {
namespace fold {

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

namespace {

using internal::gctd;
using internal::gp;

// Common:
std::vector<std::pair<int, int>> base_pairs;
std::vector<int> branch_count;

// MFE and suboptimal folding:
int max_structures;
std::multiset<computed_t, computed_energy_comparator_t> best_computeds;

// Partition function:
const int PT_MAX_BITS = 6;
const int CTD_MAX_BITS = 4;
constexpr int PT_MASK = (1 << PT_MAX_BITS) - 1;
constexpr int CTD_MASK = (1 << CTD_MAX_BITS) - 1;

struct substructure_id_t {
  constexpr static int BITS = (PT_MAX_BITS + CTD_MAX_BITS) * (1 << PT_MAX_BITS);
  constexpr static int BYTES = BITS / 8 + (BITS % 8 ? 1 : 0);
  uint16_t bits[BYTES / 2];

  bool operator<(const substructure_id_t& o) const {
    return memcmp(&bits, &o.bits, std::size_t(BYTES)) < 0;
  }
};

bool compute_partition;
partition::partition_t partition;
partition::probabilities_t probabilities;
struct nothing_t {};
SplayMap<substructure_id_t, nothing_t> substructure_map;

substructure_id_t WriteBits(int st, int en, int N, bool inside) {
  static_assert(PT_MAX_BITS + CTD_MAX_BITS <= 16, "substructure block does not fit in uint16_t");
  static_assert((-1 & PT_MASK) == PT_MASK, "mfw not a two's complement machine");
  substructure_id_t struc = {};  // Zero initialise.
  for (int i = 0, b = 0; i < N; ++i, b += PT_MAX_BITS + CTD_MAX_BITS) {
    if (inside && (i < st || i > en)) continue;
    if (!inside && i > st && i < en) continue;
    uint16_t pack = uint16_t((gp[i] & PT_MASK) << CTD_MAX_BITS | (gctd[i] & CTD_MASK));

    int byte = b / 16;
    int bit = b & 15;
    struc.bits[byte] = uint16_t(struc.bits[byte] | (pack << bit));
    int space = 16 - bit;
    if (space < CTD_MAX_BITS + PT_MAX_BITS)
      struc.bits[byte + 1] = uint16_t(struc.bits[byte + 1] | (pack >> space));
  }
  return struc;
}

substructure_id_t BuildInsideStructure(int st, int en, int N) {
  // Don't include the ctd value at st, since that's for the outside.
  return WriteBits(st + 1, en, N, true);
}

substructure_id_t BuildOutsideStructure(int st, int en, int N) {
  // Don't include the ctd value at en, since that's for the inside.
  return WriteBits(st, en + 1, N, false);
}

void AddAllCombinations(int idx) {
  const int N = static_cast<int>(gr.size());
  // Base case
  if (idx == N) {
    auto computed = energy::ComputeEnergyWithCtds({{gr, gp}, gctd, 0}, energy::gem);
    if (compute_partition) {
      penergy_t boltzmann = energy::Boltzmann(computed.energy);
      partition.q += boltzmann;
      for (int i = 0; i < N; ++i) {
        if (i < gp[i]) {
          const auto inside_structure = BuildInsideStructure(i, gp[i], N);
          const auto outside_structure = BuildOutsideStructure(i, gp[i], N);
          const bool inside_new = !substructure_map.Find(inside_structure);
          const bool outside_new = !substructure_map.Find(outside_structure);
          if (inside_new || outside_new) {
            energy_t inside_energy = energy::ComputeSubstructureEnergy(
                computed, false, i, gp[i], energy::gem);  // TODO optimisation?
            if (inside_new) {
              partition.p[i][gp[i]][0] += energy::Boltzmann(inside_energy);
              substructure_map.Insert(inside_structure, nothing_t());
            }
            if (outside_new) {
              partition.p[gp[i]][i][0] += energy::Boltzmann(computed.energy - inside_energy);
              substructure_map.Insert(outside_structure, nothing_t());
            }
          }
          probabilities[i][gp[i]][0] += boltzmann;
        }
      }
    } else {
      if (static_cast<int>(best_computeds.size()) < max_structures ||
          best_computeds.rbegin()->energy > computed.energy)
        best_computeds.insert(std::move(computed));
      if (static_cast<int>(best_computeds.size()) > max_structures)
        best_computeds.erase(--best_computeds.end());
    }
    return;
  }

  // If we already set this, this isn't a valid base pair, it's not part of a multiloop, can't set
  // ctds so continue.
  if (gctd[idx] != CTD_NA || gp[idx] == -1 || branch_count[idx] < 2) {
    AddAllCombinations(idx + 1);
    return;
  }

  const bool lu_exists = idx - 1 >= 0 && gp[idx - 1] == -1;
  const bool lu_shared = lu_exists && idx - 2 >= 0 && gp[idx - 2] != -1;
  const bool lu_usable = lu_exists &&
      (!lu_shared ||
          (gctd[gp[idx - 2]] != CTD_3_DANGLE && gctd[gp[idx - 2]] != CTD_MISMATCH &&
              gctd[gp[idx - 2]] != CTD_RCOAX_WITH_PREV));
  const bool ru_exists = gp[idx] + 1 < N && gp[gp[idx] + 1] == -1;
  const bool ru_shared = ru_exists && gp[idx] + 2 < N && gp[gp[idx] + 2] != -1;
  const bool ru_usable = ru_exists &&
      (!ru_shared ||
          (gctd[gp[idx] + 2] != CTD_5_DANGLE && gctd[gp[idx] + 2] != CTD_MISMATCH &&
              gctd[gp[idx] + 2] != CTD_LCOAX_WITH_NEXT));
  // Even if the next branch is an outer branch, everything will be magically handled.
  // CTD_UNUSED
  gctd[idx] = CTD_UNUSED;
  AddAllCombinations(idx + 1);

  // CTD_3_DANGLE
  if (ru_usable) {
    gctd[idx] = CTD_3_DANGLE;
    AddAllCombinations(idx + 1);
  }

  // CTD_5_DANGLE
  if (lu_usable) {
    gctd[idx] = CTD_5_DANGLE;
    AddAllCombinations(idx + 1);
  }

  // CTD_MISMATCH
  if (ru_usable && lu_usable) {
    gctd[idx] = CTD_MISMATCH;
    AddAllCombinations(idx + 1);
  }

  // Check that the next branch hasn't been set already. If it's unused or na, try re-writing it.
  // CTD_LCOAX_WITH_NEXT
  if (lu_usable && ru_usable && ru_shared) {
    auto prevval = gctd[gp[idx] + 2];
    if (prevval == CTD_UNUSED || prevval == CTD_NA) {
      gctd[idx] = CTD_LCOAX_WITH_NEXT;
      gctd[gp[idx] + 2] = CTD_LCOAX_WITH_PREV;
      AddAllCombinations(idx + 1);
      gctd[gp[idx] + 2] = prevval;
    }
  }

  // Check that the previous branch hasn't been set already.
  // CTD_RCOAX_WITH_PREV
  if (lu_usable && lu_shared && ru_usable) {
    auto prevval = gctd[gp[idx - 2]];
    if (prevval == CTD_UNUSED || prevval == CTD_NA) {
      gctd[idx] = CTD_RCOAX_WITH_PREV;
      gctd[gp[idx - 2]] = CTD_RCOAX_WITH_NEXT;
      AddAllCombinations(idx + 1);
      gctd[gp[idx - 2]] = prevval;
    }
  }

  // CTD_FCOAX_WITH_NEXT
  if (gp[idx] + 1 < N && gp[gp[idx] + 1] != -1) {
    auto prevval = gctd[gp[idx] + 1];
    if (prevval == CTD_UNUSED || prevval == CTD_NA) {
      gctd[idx] = CTD_FCOAX_WITH_NEXT;
      gctd[gp[idx] + 1] = CTD_FCOAX_WITH_PREV;
      AddAllCombinations(idx + 1);
      gctd[gp[idx] + 1] = prevval;
    }
  }

  // Reset back to NA.
  gctd[idx] = CTD_NA;
}

void BruteForce(int idx) {
  if (idx == static_cast<int>(base_pairs.size())) {
    // Small optimisation for case when we're just getting one structure.
    if (max_structures == 1 && !compute_partition) {
      auto computed = energy::ComputeEnergy({gr, gp}, energy::gem);
      if (best_computeds.empty() || computed.energy < best_computeds.begin()->energy)
        best_computeds.insert(std::move(computed));
      if (best_computeds.size() == 2) best_computeds.erase(--best_computeds.end());
    } else {
      // Precompute whether things are multiloops or not.
      branch_count = internal::GetBranchCounts(gp);
      AddAllCombinations(0);
    }

    return;
  }
  // Don't take this base pair.
  BruteForce(idx + 1);

  // Take this base pair.
  bool can_take = true;
  const auto& pair = base_pairs[idx];
  // Only need to check in the range of this base pair. Since we ordered by
  // increasing st, anything at or after this will either be the start of something starting at st,
  // or something ending, both of which conflict with this base pair.
  for (int i = pair.first; i <= pair.second; ++i) {
    if (gp[i] != -1) {
      can_take = false;
      break;
    }
  }
  if (can_take) {
    gp[pair.first] = pair.second;
    gp[pair.second] = pair.first;
    BruteForce(idx + 1);
    gp[pair.first] = -1;
    gp[pair.second] = -1;
  }
}

void InvokeBruteForce(const primary_t& r, const energy::EnergyModel& em, int max_structures_,
    bool compute_partition_, bool allow_lonely_pairs) {
  SetFoldGlobalState(r, em);
  best_computeds.clear();
  base_pairs.clear();
  max_structures = max_structures_ == -1 ? internal::MAX_STRUCTURES : max_structures_;
  compute_partition = compute_partition_;
  if (compute_partition_) {
    partition.q = 0;
    partition.p = array3d_t<penergy_t, 1>(r.size(), 0);
    probabilities = array3d_t<penergy_t, 1>(r.size(), 0);
  }
  // Add base pairs in order of increasing st, then en.
  for (int st = 0; st < static_cast<int>(r.size()); ++st) {
    for (int en = st + HAIRPIN_MIN_SZ + 1; en < static_cast<int>(r.size()); ++en) {
      bool allowed = allow_lonely_pairs ? CanPair(r[st], r[en]) : energy::ViableFoldingPair(st, en);
      if (allowed) base_pairs.emplace_back(st, en);
    }
  }
  BruteForce(0);
}

}  // namespace

computed_t FoldBruteForce(const primary_t& r, const energy::EnergyModel& em) {
  return SuboptimalBruteForce(r, em, 1)[0];
}

std::vector<computed_t> SuboptimalBruteForce(
    const primary_t& r, const energy::EnergyModel& em, int max_structures_) {
  InvokeBruteForce(r, em, max_structures_, false, false);
  return std::vector<computed_t>(best_computeds.begin(), best_computeds.end());
}

std::pair<partition::partition_t, partition::probabilities_t> PartitionBruteForce(
    const primary_t& r, const energy::EnergyModel& em) {
  const int N = static_cast<int>(r.size());
  // Preconditions:
  static_assert(CTD_SIZE < (1 << CTD_MAX_BITS), "need increase ctd bits for brute force");
  // Plus one to N, since -1 takes up a spot.
  verify(N + 1 < (1 << PT_MAX_BITS), "sequence too long for brute force partition");

  InvokeBruteForce(r, em, 1, true, false);  // Allow lonely pairs for the partition function. TODO?
  substructure_map.Clear();  // Don't waste memory.
  for (int i = 0; i < N; ++i)
    for (int j = i; j < N; ++j) probabilities[i][j][0] /= partition.q;
  return {std::move(partition), std::move(probabilities)};
}

}  // namespace fold
}  // namespace memerna
