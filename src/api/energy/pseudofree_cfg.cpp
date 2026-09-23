// Copyright 2025 Eliot Courtney.
#include "api/energy/pseudofree_cfg.h"

#include <algorithm>
#include <cstddef>
#include <ostream>
#include <span>
#include <utility>
#include <vector>

#include "model/energy.h"
#include "model/primary.h"
#include "util/argparse.h"
#include "util/error.h"

namespace mrna::erg {

namespace {

std::vector<Energy> BuildUnpairedSum(const std::vector<Energy>& unpaired) {
  if (unpaired.empty()) return {};
  std::vector<Energy> unpaired_sum;
  unpaired_sum.reserve(unpaired.size() + 1);
  unpaired_sum.push_back(ZERO_E);
  for (const auto& e : unpaired) unpaired_sum.push_back(unpaired_sum.back() + e);
  return unpaired_sum;
}

std::vector<BoltzEnergy> BuildBoltz(std::span<const Energy> energies) {
  std::vector<BoltzEnergy> boltz;
  boltz.resize(energies.size());
  for (size_t i = 0; i < energies.size(); ++i) boltz[i] = Boltz(energies[i]);
  return boltz;
}

std::vector<BoltzEnergy> BuildUnpairedSumLog(std::span<const Energy> unpaired) {
  if (unpaired.empty()) return {};
  std::vector<BoltzEnergy> unpaired_sum_log;
  unpaired_sum_log.reserve(unpaired.size() + 1);
  unpaired_sum_log.push_back(ZERO_B);
  for (const auto& e : unpaired) {
    const auto log_boltz = e >= CAP_E ? ZERO_B : e.LogBoltz();
    unpaired_sum_log.push_back(unpaired_sum_log.back() + log_boltz);
  }
  return unpaired_sum_log;
}

std::vector<std::size_t> BuildUnpairedZeroCount(std::span<const Energy> unpaired) {
  if (unpaired.empty()) return {};
  std::vector<std::size_t> zero_count;
  zero_count.reserve(unpaired.size() + 1);
  zero_count.push_back(0);
  for (const auto& e : unpaired)
    zero_count.push_back(zero_count.back() + static_cast<std::size_t>(e >= CAP_E));
  return zero_count;
}

}  // namespace

void RegisterOptsPseudofree(ArgParse* args) {
  args->RegisterOpt(OPT_PAIRED_PSEUDOFREE);
  args->RegisterOpt(OPT_UNPAIRED_PSEUDOFREE);
}

PseudofreeCfg::PseudofreeCfg(std::vector<Energy> paired, std::vector<Energy> unpaired)
    : paired_(std::move(paired)), unpaired_(std::move(unpaired)),
      unpaired_sum_(BuildUnpairedSum(unpaired_)) {}

Energy PseudofreeCfg::MinUnpaired() const {
  if (unpaired_.empty()) return ZERO_E;
  Energy min_unpaired = unpaired_[0];
  for (const auto& e : unpaired_) min_unpaired = std::min(min_unpaired, e);
  return min_unpaired;
}

Energy PseudofreeCfg::SmallestSubsetSum() const {
  Energy sum = ZERO_E;
  for (const auto& e : unpaired_)
    if (e < ZERO_E) sum += e;
  for (const auto& e : paired_)
    if (e < ZERO_E) sum += e;
  return sum;
}

void PseudofreeCfg::Verify(const Primary& r) const {
  if (!paired_.empty())
    verify(paired_.size() == std::size_t(r.size()), "pseudofree paired must be same length as seq");
  if (!unpaired_.empty())
    verify(unpaired_.size() == std::size_t(r.size()),
        "pseudofree unpaired must be same length as seq");
}

PseudofreeCfg PseudofreeCfg::FromArgParse(const ArgParse& args) {
  auto pf_paired = args.GetMultipleOr<Energy>(OPT_PAIRED_PSEUDOFREE, {});
  auto pf_unpaired = args.GetMultipleOr<Energy>(OPT_UNPAIRED_PSEUDOFREE, {});
  return {std::move(pf_paired), std::move(pf_unpaired)};
}

BoltzPseudofreeCfg PseudofreeCfg::Boltz() const { return BoltzPseudofreeCfg(*this); }

BoltzPseudofreeCfg::BoltzPseudofreeCfg(const PseudofreeCfg& pf)
    : paired_(BuildBoltz(pf.PairedEnergies())), unpaired_(BuildBoltz(pf.UnpairedEnergies())),
      unpaired_sum_log_(BuildUnpairedSumLog(pf.UnpairedEnergies())),
      unpaired_zero_count_(BuildUnpairedZeroCount(pf.UnpairedEnergies())) {}

void BoltzPseudofreeCfg::Verify(const Primary& r) const {
  if (!paired_.empty())
    verify(paired_.size() == std::size_t(r.size()), "pseudofree paired must be same length as seq");
  if (!unpaired_.empty())
    verify(unpaired_.size() == std::size_t(r.size()),
        "pseudofree unpaired must be same length as seq");
}

std::ostream& operator<<(std::ostream& str, const PseudofreeCfg& o) {
  return str << "PseudofreeCfg{paired=" << o.PairedEnergies().size()
             << ", unpaired=" << o.UnpairedEnergies().size() << "}";
}

}  // namespace mrna::erg
