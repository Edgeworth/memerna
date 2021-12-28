// Copyright 2016 E.
#ifndef COMPUTE_BRUTE_H_
#define COMPUTE_BRUTE_H_

#include <set>
#include <stack>
#include <utility>
#include <vector>

#include "compute/energy/globals.h"
#include "compute/energy/model.h"
#include "compute/energy/structure.h"
#include "compute/mfe/globals.h"
#include "compute/mfe/mfe.h"
#include "compute/partition/partition.h"
#include "model/parsing.h"
#include "model/structure.h"
#include "util/splaymap.h"

namespace mrna {

namespace internal {

std::vector<int> GetBranchCounts(const std::vector<int>& p);

}

class BruteForce {
 public:
  struct Result {
    // Common:
    std::vector<std::pair<int, int>> base_pairs;
    std::vector<int> branch_count;

    // MFE and suboptimal folding:
    int max_structures;
    std::multiset<Computed, ComputedEnergyCmp> best_computeds;

    // TODO: Switch to optional?
    bool compute_partition;
    partition::Partition partition;
    partition::Probabilities probabilities;
  };

  Result Run(const Primary& r, const energy::EnergyModel& em, int max_structures,
      bool compute_partition, bool allow_lonely_pairs);

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

    bool operator<(const SubstructureId& o) const {
      return memcmp(&bits, &o.bits, std::size_t(BYTES)) < 0;
    }
  };

  Result res_;
  SplaySet<SubstructureId> substructure_map;

  void AddAllCombinations(int idx);
  void Dfs(int idx);

  static SubstructureId WriteBits(int st, int en, int N, bool inside);
  static SubstructureId BuildInsideStructure(int st, int en, int N);
  static SubstructureId BuildOutsideStructure(int st, int en, int N);
};

}  // namespace mrna

#endif  // COMPUTE_BRUTE_H_
