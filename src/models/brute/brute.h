// Copyright 2016 Eliot Courtney.
#ifndef COMPUTE_BRUTE_BRUTE_H_
#define COMPUTE_BRUTE_BRUTE_H_

#include <compare>
#include <cstdint>
#include <set>
#include <utility>
#include <vector>

#include "api/brute/brute_cfg.h"
#include "api/energy/model.h"
#include "api/part.h"
#include "api/subopt/subopt.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "models/t04/part/part.h"
#include "util/splaymap.h"

namespace mrna::md::brute {

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
  brute::BruteCfg cfg_;

  Secondary s_;
  Ctds ctd_;
  brute::BruteResult res_{};
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

#endif  // COMPUTE_BRUTE_BRUTE_H_
