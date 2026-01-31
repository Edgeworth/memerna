// Copyright 2025 Eliot Courtney.
#ifndef API_ENERGY_PSEUDOFREE_H_
#define API_ENERGY_PSEUDOFREE_H_

#include <vector>

#include "model/energy.h"
#include "model/primary.h"
#include "util/argparse.h"

namespace mrna::erg {

inline const Opt OPT_PAIRED_PSEUDOFREE =
    Opt(Opt::ARG)
        .LongName("pf-paired")
        .Multiple()
        .Help("comma separated energies for paired pseudofree energy");
inline const Opt OPT_UNPAIRED_PSEUDOFREE =
    Opt(Opt::ARG)
        .LongName("pf-unpaired")
        .Multiple()
        .Help("comma separated energies for unpaired pseudofree energy");

void RegisterOptsPseudofree(ArgParse* args);

// Pseudofree energy configuration with precomputed cumulative sums.
// Pseudofree energies are per-position adjustments to the energy model.
class PseudofreeCfg {
 public:
  // Pseudofree energies. Ignored if empty.
  std::vector<Energy> paired;
  std::vector<Energy> unpaired;
  // Cumulative sum of size N+1 (first element is nothing).
  std::vector<Energy> unpaired_cum;

  [[nodiscard]] constexpr Energy Unpaired(int n) const {
    if (unpaired.empty()) return ZERO_E;
    return unpaired[n];
  }

  // Inclusive range, unlike pf_unpaired_cum directly.
  [[nodiscard]] constexpr Energy UnpairedCum(int st, int en) const {
    assert(st <= en + 1);
    if (unpaired.empty()) return ZERO_E;
    return unpaired_cum[en + 1] - unpaired_cum[st];
  }

  [[nodiscard]] constexpr Energy Paired(int st, int en) const {
    assert(st <= en + 1);
    if (paired.empty()) return ZERO_E;
    return paired[st] + paired[en];
  }

  void Load(std::vector<Energy> pf_paired, std::vector<Energy> pf_unpaired);
  void Verify(const Primary& r) const;

  static PseudofreeCfg FromArgParse(const ArgParse& args);
};

class BoltzPseudofreeCfg {
 public:
  // Pseudofree energies. Ignored if empty.
  std::vector<BoltzEnergy> paired;
  std::vector<BoltzEnergy> unpaired;
  // Cumulative sum of size N+1 (first element is nothing).
  std::vector<BoltzEnergy> unpaired_cum_log;

  [[nodiscard]] BoltzEnergy Unpaired(int n) const {
    if (unpaired.empty()) return ONE_B;
    return unpaired[n];
  }

  // Inclusive range, unlike pf_unpaired_cum directly.
  [[nodiscard]] BoltzEnergy UnpairedCum(int st, int en) const {
    assert(st <= en + 1);
    if (unpaired.empty()) return ONE_B;
    return exp(unpaired_cum_log[en + 1] - unpaired_cum_log[st]);
  }

  [[nodiscard]] BoltzEnergy Paired(int st, int en) const {
    assert(st <= en + 1);
    if (paired.empty()) return ONE_B;
    return paired[st] * paired[en];
  }

  void Load(std::vector<Energy> pf_paired, std::vector<Energy> pf_unpaired);
  void Verify(const Primary& r) const;
};

}  // namespace mrna::erg

#endif  // API_ENERGY_PSEUDOFREE_H_
