// Copyright 2016 Eliot Courtney.
#ifndef BACKENDS_BRUTE_BRUTE_H_
#define BACKENDS_BRUTE_BRUTE_H_

#include <cstdint>
#include <set>
#include <utility>
#include <vector>

#include "api/brute/brute_cfg.h"
#include "api/ctx/backend.h"
#include "api/energy/energy_cfg.h"
#include "api/subopt/subopt.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/pfn.h"
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
  // Subopt structures:
  // TODO(2): use splayset here?
  std::multiset<subopt::SuboptResult, SuboptCmp> subopts;

  // Pfn function:
  PfnTables pfn;
};

class Brute {
 public:
  Brute(const Primary& r, BackendModelPtr m, BruteCfg cfg);

  BruteResult Run();

 private:
  // Pfn function:
  constexpr static uint32_t PT_MAX_BITS = 5;
  constexpr static uint32_t CTD_MAX_BITS = 4;
  constexpr static uint32_t PT_MASK = (1 << PT_MAX_BITS) - 1;
  constexpr static uint32_t CTD_MASK = (1 << CTD_MAX_BITS) - 1;

  struct SubstructureId {
    constexpr static uint32_t BITS = (PT_MAX_BITS + CTD_MAX_BITS) * (1 << PT_MAX_BITS);
    constexpr static uint32_t BYTES = BITS / 8 + (BITS % 8 ? 1 : 0);
    uint16_t bits[BYTES / 2];

    constexpr auto operator<=>(const SubstructureId&) const = default;
  };

  Primary r_;
  BackendModelPtr m_;
  BackendBoltzModelPtr bm_;
  BackendModelPtr underlying_;
  erg::EnergyCfg energy_cfg_;
  BruteCfg brute_cfg_;

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

#endif  // BACKENDS_BRUTE_BRUTE_H_
