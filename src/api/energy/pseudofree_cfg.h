// Copyright 2025 Eliot Courtney.
#ifndef API_ENERGY_PSEUDOFREE_CFG_H_
#define API_ENERGY_PSEUDOFREE_CFG_H_

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <iosfwd>
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
  const std::vector<Energy> paired{};
  const std::vector<Energy> unpaired{};
  // Cumulative sum of size N+1 (first element is nothing).
  const std::vector<Energy> unpaired_sum{};

  PseudofreeCfg() = default;
  PseudofreeCfg(std::vector<Energy> paired, std::vector<Energy> unpaired);

  constexpr auto operator<=>(const PseudofreeCfg&) const = default;

  [[nodiscard]] constexpr bool Empty() const { return paired.empty() && unpaired.empty(); }

  [[nodiscard]] constexpr Energy Unpaired(int n) const {
    if (unpaired.empty()) return ZERO_E;
    return unpaired[n];
  }

  // Inclusive range, unlike unpaired_sum directly.
  [[nodiscard]] constexpr Energy UnpairedSum(int st, int en) const {
    assert(st <= en + 1);
    if (unpaired.empty()) return ZERO_E;
    return unpaired_sum[en + 1] - unpaired_sum[st];
  }

  [[nodiscard]] constexpr Energy Paired(int st, int en) const {
    assert(st <= en + 1);
    if (paired.empty()) return ZERO_E;
    return paired[st] + paired[en];
  }

  void Verify(const Primary& r) const;

  static PseudofreeCfg FromArgParse(const ArgParse& args);
};

std::ostream& operator<<(std::ostream& str, const PseudofreeCfg& o);

class BoltzPseudofreeCfg {
 public:
  const std::vector<BoltzEnergy> paired{};
  const std::vector<BoltzEnergy> unpaired{};
  const std::vector<BoltzEnergy> unpaired_sum_log{};

  BoltzPseudofreeCfg() = default;
  explicit BoltzPseudofreeCfg(const PseudofreeCfg& pf);

  [[nodiscard]] bool Empty() const { return paired.empty() && unpaired.empty(); }

  [[nodiscard]] BoltzEnergy Unpaired(int n) const {
    if (unpaired.empty()) return ONE_B;
    return unpaired[n];
  }

  // Inclusive range, unlike unpaired_sum_log directly.
  [[nodiscard]] BoltzEnergy UnpairedProd(int st, int en) const {
    assert(st <= en + 1);
    if (unpaired.empty()) return ONE_B;
    return exp(unpaired_sum_log[en + 1] - unpaired_sum_log[st]);
  }

  [[nodiscard]] BoltzEnergy Paired(int st, int en) const {
    assert(st <= en + 1);
    if (paired.empty()) return ONE_B;
    return paired[st] * paired[en];
  }

  void Verify(const Primary& r) const;
};

}  // namespace mrna::erg

template <>
struct fmt::formatter<mrna::erg::PseudofreeCfg> : ostream_formatter {};

#endif  // API_ENERGY_PSEUDOFREE_CFG_H_
