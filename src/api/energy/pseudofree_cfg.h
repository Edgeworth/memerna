// Copyright 2025 Eliot Courtney.
#ifndef API_ENERGY_PSEUDOFREE_CFG_H_
#define API_ENERGY_PSEUDOFREE_CFG_H_

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <iosfwd>
#include <ostream>
#include <span>
#include <vector>

#include "model/energy.h"
#include "model/primary.h"
#include "util/argparse.h"

namespace mrna::erg {

class BoltzPseudofreeCfg;
class NoBoltzPseudofreeCfg;

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
  PseudofreeCfg() = default;
  PseudofreeCfg(std::vector<Energy> paired, std::vector<Energy> unpaired);

  constexpr auto operator<=>(const PseudofreeCfg&) const = default;

  [[nodiscard]] constexpr bool Empty() const { return paired_.empty() && unpaired_.empty(); }
  [[nodiscard]] constexpr bool HasPaired() const { return !paired_.empty(); }
  [[nodiscard]] constexpr bool HasUnpaired() const { return !unpaired_.empty(); }
  [[nodiscard]] constexpr std::span<const Energy> PairedEnergies() const { return paired_; }
  [[nodiscard]] constexpr std::span<const Energy> UnpairedEnergies() const { return unpaired_; }
  [[nodiscard]] Energy MinUnpaired() const;
  [[nodiscard]] Energy SmallestSubsetSum() const;

  [[nodiscard]] constexpr Energy Unpaired(int n) const {
    if (unpaired_.empty()) return ZERO_E;
    assert(n >= 0 && n < static_cast<int>(unpaired_.size()));
    return unpaired_[n];
  }

  // Inclusive range over unpaired positions.
  [[nodiscard]] constexpr Energy UnpairedSum(int st, int en) const {
    assert(st <= en + 1);
    if (unpaired_.empty()) return ZERO_E;
    assert(st >= 0 && en + 1 < static_cast<int>(unpaired_sum_.size()));
    return unpaired_sum_[en + 1] - unpaired_sum_[st];
  }

  [[nodiscard]] constexpr Energy Paired(int st, int en) const {
    assert(st <= en);
    if (paired_.empty()) return ZERO_E;
    assert(st >= 0 && en < static_cast<int>(paired_.size()));
    return paired_[st] + paired_[en];
  }

  void Verify(const Primary& r) const;

  [[nodiscard]] BoltzPseudofreeCfg Boltz() const;

  static PseudofreeCfg FromArgParse(const ArgParse& args);

 private:
  // Pseudofree energies. Ignored if empty.
  const std::vector<Energy> paired_{};
  const std::vector<Energy> unpaired_{};
  // Cumulative sum of size N+1 (first element is nothing).
  const std::vector<Energy> unpaired_sum_{};
};

class NoPseudofreeCfg {
 public:
  constexpr NoPseudofreeCfg() = default;
  constexpr explicit NoPseudofreeCfg([[maybe_unused]] const PseudofreeCfg& pf) {}

  constexpr auto operator<=>(const NoPseudofreeCfg&) const = default;

  [[nodiscard]] static constexpr bool Empty() { return true; }
  [[nodiscard]] static constexpr bool HasPaired() { return false; }
  [[nodiscard]] static constexpr bool HasUnpaired() { return false; }
  [[nodiscard]] static constexpr std::span<const Energy> PairedEnergies() { return {}; }
  [[nodiscard]] static constexpr std::span<const Energy> UnpairedEnergies() { return {}; }
  [[nodiscard]] static constexpr Energy MinUnpaired() { return ZERO_E; }
  [[nodiscard]] static constexpr Energy SmallestSubsetSum() { return ZERO_E; }
  [[nodiscard]] static constexpr Energy Unpaired([[maybe_unused]] int n) { return ZERO_E; }
  [[nodiscard]] static constexpr Energy UnpairedSum(
      [[maybe_unused]] int st, [[maybe_unused]] int en) {
    return ZERO_E;
  }
  [[nodiscard]] static constexpr Energy Paired([[maybe_unused]] int st, [[maybe_unused]] int en) {
    return ZERO_E;
  }

  static void Verify([[maybe_unused]] const Primary& r) {}

  [[nodiscard]] static constexpr NoBoltzPseudofreeCfg Boltz();
};

std::ostream& operator<<(std::ostream& str, const PseudofreeCfg& o);

class BoltzPseudofreeCfg {
 public:
  BoltzPseudofreeCfg() = default;
  explicit BoltzPseudofreeCfg(const PseudofreeCfg& pf);
  explicit BoltzPseudofreeCfg([[maybe_unused]] const NoPseudofreeCfg& pf) {}

  [[nodiscard]] const BoltzPseudofreeCfg& Boltz() const { return *this; }

  [[nodiscard]] bool Empty() const { return paired_.empty() && unpaired_.empty(); }

  [[nodiscard]] BoltzEnergy Unpaired(int n) const {
    if (unpaired_.empty()) return ONE_B;
    assert(n >= 0 && n < static_cast<int>(unpaired_.size()));
    return unpaired_[n];
  }

  // Inclusive range over unpaired positions.
  [[nodiscard]] BoltzEnergy UnpairedProd(int st, int en) const {
    assert(st <= en + 1);
    if (unpaired_.empty()) return ONE_B;
    assert(st >= 0 && en + 1 < static_cast<int>(unpaired_sum_log_.size()));
    return exp(unpaired_sum_log_[en + 1] - unpaired_sum_log_[st]);
  }

  [[nodiscard]] BoltzEnergy Paired(int st, int en) const {
    assert(st <= en);
    if (paired_.empty()) return ONE_B;
    assert(st >= 0 && en < static_cast<int>(paired_.size()));
    return paired_[st] * paired_[en];
  }

  void Verify(const Primary& r) const;

 private:
  const std::vector<BoltzEnergy> paired_{};
  const std::vector<BoltzEnergy> unpaired_{};
  const std::vector<BoltzEnergy> unpaired_sum_log_{};
};

class NoBoltzPseudofreeCfg {
 public:
  constexpr NoBoltzPseudofreeCfg() = default;
  constexpr explicit NoBoltzPseudofreeCfg([[maybe_unused]] const NoPseudofreeCfg& pf) {}

  [[nodiscard]] static constexpr NoBoltzPseudofreeCfg Boltz() { return {}; }

  [[nodiscard]] static constexpr bool Empty() { return true; }
  [[nodiscard]] static BoltzEnergy Unpaired([[maybe_unused]] int n) { return ONE_B; }
  [[nodiscard]] static BoltzEnergy UnpairedProd([[maybe_unused]] int st, [[maybe_unused]] int en) {
    return ONE_B;
  }
  [[nodiscard]] static BoltzEnergy Paired([[maybe_unused]] int st, [[maybe_unused]] int en) {
    return ONE_B;
  }

  static void Verify([[maybe_unused]] const Primary& r) {}
};

[[nodiscard]] constexpr NoBoltzPseudofreeCfg NoPseudofreeCfg::Boltz() { return {}; }

inline std::ostream& operator<<(std::ostream& str, [[maybe_unused]] const NoPseudofreeCfg& o) {
  return str << "{}";
}

}  // namespace mrna::erg

template <>
struct fmt::formatter<mrna::erg::PseudofreeCfg> : ostream_formatter {};

template <>
struct fmt::formatter<mrna::erg::NoPseudofreeCfg> : ostream_formatter {};

#endif  // API_ENERGY_PSEUDOFREE_CFG_H_
