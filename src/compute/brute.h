// Copyright 2016 E.
#ifndef COMPUTE_BRUTE_H_
#define COMPUTE_BRUTE_H_

#include <cstdint>
#include <set>
#include <utility>
#include <vector>

#include "compute/boltz_dp.h"
#include "compute/energy/model.h"
#include "compute/partition/partition.h"
#include "compute/subopt/subopt.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/splaymap.h"

namespace mrna {

class BruteForce {
 public:
  struct Result {
   private:
    struct SuboptCmp {
      bool operator()(const subopt::SuboptResult& a, const subopt::SuboptResult& b) const {
        // Kept in a multiset, so this is just used for ordering, not deduplication.
        // There should be no duplicates added anyway. Single comparison to keep it fast.
        return a.energy < b.energy;
      }
    };

   public:
    // Common:
    std::vector<std::pair<int, int>> base_pairs;
    std::vector<int> branch_count;

    // MFE and suboptimal folding:
    int max_structures;
    std::multiset<subopt::SuboptResult, SuboptCmp> subopts;

    // TODO: Switch to optional?
    bool compute_partition;
    partition::Partition partition;
    Probabilities probabilities;
  };

  Result Run(Primary r, const energy::EnergyModel& em, int max_structures, bool compute_partition,
      bool allow_lonely_pairs);

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
  Secondary s_;
  Ctds ctd_;
  energy::EnergyModel em_;
  Result res_;
  SplaySet<SubstructureId> substructure_map_;

  void AddAllCombinations(int idx);
  void Dfs(int idx);

  SubstructureId WriteBits(int st, int en, int N, bool inside);
  SubstructureId BuildInsideStructure(int st, int en, int N);
  SubstructureId BuildOutsideStructure(int st, int en, int N);
};

}  // namespace mrna

#endif  // COMPUTE_BRUTE_H_
