// Copyright 2016 E.
#ifndef MODELS_BRUTE_BRUTE_H_
#define MODELS_BRUTE_BRUTE_H_

#include <compare>
#include <cstdint>
#include <set>
#include <utility>
#include <vector>

#include "api/brute/brute_cfg.h"
#include "api/energy/energy_cfg.h"
#include "api/energy/model.h"
#include "api/subopt/subopt.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/part.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/splaymap.h"

namespace mrna::md::brute {

using mrna::brute::BruteCfg;

struct SuboptCmp {
  bool operator()(const subopt::SuboptResult& a, const subopt::SuboptResult& b) const {
    // Kept in a multiset, so this is just used for ordering, not deduplication.
    // There should be no duplicates added anyway. Single comparison to keep it fast.
    return a.energy < b.energy;
  }
};

struct BruteResult {
  // Suboptimal structures:
  // TODO(2): use splayset here?
  std::multiset<subopt::SuboptResult, SuboptCmp> subopts;

  // Partition function:
  Part part;
};

class Brute {
 public:
  Brute(const Primary& r, erg::EnergyModelPtr em, BruteCfg cfg);

  BruteResult Run();

 private:
  // Partition function:
  inline constexpr static uint32_t PT_MAX_BITS = 6;
  inline constexpr static uint32_t CTD_MAX_BITS = 4;
  inline constexpr static uint32_t PT_MASK = (1 << PT_MAX_BITS) - 1;
  inline constexpr static uint32_t CTD_MASK = (1 << CTD_MAX_BITS) - 1;

  struct SubstructureId {
    inline constexpr static uint32_t BITS = (PT_MAX_BITS + CTD_MAX_BITS) * (1 << PT_MAX_BITS);
    inline constexpr static uint32_t BYTES = BITS / 8 + (BITS % 8 ? 1 : 0);
    uint16_t bits[BYTES / 2];

    constexpr auto operator<=>(const SubstructureId&) const = default;
  };

  Primary r_;
  erg::EnergyModelPtr em_;
  erg::BoltzEnergyModelPtr bem_;
  erg::EnergyModelPtr underlying_;
  erg::EnergyCfg::Ctd ctd_cfg_;
  BruteCfg cfg_;

  Secondary s_;
  Ctds ctd_;
  BruteResult res_{};
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

}  // namespace mrna::md::brute

#endif  // MODELS_BRUTE_BRUTE_H_
