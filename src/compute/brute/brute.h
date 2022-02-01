// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_BRUTE_BRUTE_H_
#define COMPUTE_BRUTE_BRUTE_H_

#include <compare>
#include <cstdint>
#include <set>
#include <utility>
#include <vector>

#include "compute/boltz_dp.h"
#include "compute/brute/config.h"
#include "compute/energy/model.h"
#include "compute/partition/partition.h"
#include "compute/subopt/subopt.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/splaymap.h"

namespace mrna::brute {

struct SuboptCmp {
  bool operator()(const subopt::SuboptResult& a, const subopt::SuboptResult& b) const {
    // Kept in a multiset, so this is just used for ordering, not deduplication.
    // There should be no duplicates added anyway. Single comparison to keep it fast.
    return a.energy < b.energy;
  }
};

struct BruteResult {
  // Suboptimal structures:
  // TODO: use splayset here?
  std::multiset<subopt::SuboptResult, SuboptCmp> subopts;

  // Partition function:
  part::Part part;
  BoltzProbs prob;
};

class BruteForce {
 public:
  BruteForce(const Primary& r, energy::EnergyModelPtr em, BruteCfg cfg);

  BruteResult Run();

 private:
  // Partition function:
  inline constexpr static int PT_MAX_BITS = 6;
  inline constexpr static int CTD_MAX_BITS = 4;
  inline constexpr static int PT_MASK = (1 << PT_MAX_BITS) - 1;
  inline constexpr static int CTD_MASK = (1 << CTD_MAX_BITS) - 1;

  struct SubstructureId {
    inline constexpr static int BITS = (PT_MAX_BITS + CTD_MAX_BITS) * (1 << PT_MAX_BITS);
    inline constexpr static int BYTES = BITS / 8 + (BITS % 8 ? 1 : 0);
    uint16_t bits[BYTES / 2];

    auto operator<=>(const SubstructureId&) const = default;
  };

  Primary r_;
  energy::EnergyModelPtr em_;
  BruteCfg cfg_;

  Secondary s_;
  Ctds ctd_;
  BruteResult res_;
  std::vector<std::pair<int, int>> pairs_;  // Holds all possible base pairs to try.
  std::vector<int> branch_count_;  // Number of sibling branches.

  SplaySet<SubstructureId> substructure_map_;

  void Dfs(int idx);
  void AddAllCombinations(int idx);

  void PruneInsertSubopt(Energy e);

  SubstructureId WriteBits(int st, int en, int N, bool inside);
  SubstructureId BuildInsideStructure(int st, int en, int N);
  SubstructureId BuildOutsideStructure(int st, int en, int N);
};

}  // namespace mrna::brute

#endif  // COMPUTE_BRUTE_BRUTE_H_
