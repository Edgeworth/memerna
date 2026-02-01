// Copyright 2025 Eliot Courtney.
#include "api/energy/pseudofree_cfg.h"

#include <utility>
#include <vector>

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

std::vector<BoltzEnergy> BuildBoltz(const std::vector<Energy>& energies) {
  std::vector<BoltzEnergy> boltz;
  boltz.resize(energies.size());
  for (size_t i = 0; i < energies.size(); ++i) boltz[i] = energies[i].Boltz();
  return boltz;
}

std::vector<BoltzEnergy> BuildUnpairedSumLog(const std::vector<Energy>& unpaired) {
  if (unpaired.empty()) return {};
  std::vector<BoltzEnergy> unpaired_sum_log;
  unpaired_sum_log.reserve(unpaired.size() + 1);
  unpaired_sum_log.push_back(ZERO_B);
  for (const auto& e : unpaired) unpaired_sum_log.push_back(unpaired_sum_log.back() + e.LogBoltz());
  return unpaired_sum_log;
}

}  // namespace

void RegisterOptsPseudofree(ArgParse* args) {
  args->RegisterOpt(OPT_PAIRED_PSEUDOFREE);
  args->RegisterOpt(OPT_UNPAIRED_PSEUDOFREE);
}

PseudofreeCfg::PseudofreeCfg(std::vector<Energy> paired_, std::vector<Energy> unpaired_)
    : paired(std::move(paired_)), unpaired(std::move(unpaired_)),
      unpaired_sum(BuildUnpairedSum(unpaired)) {}

void PseudofreeCfg::Verify(const Primary& r) const {
  if (!paired.empty())
    verify(paired.size() == r.size(), "pseudofree paired must be same length as seq");
  if (!unpaired.empty())
    verify(unpaired.size() == r.size(), "pseudofree unpaired must be same length as seq");
}

PseudofreeCfg PseudofreeCfg::FromArgParse(const ArgParse& args) {
  auto pf_paired = args.GetMultipleOr<Energy>(OPT_PAIRED_PSEUDOFREE, {});
  auto pf_unpaired = args.GetMultipleOr<Energy>(OPT_UNPAIRED_PSEUDOFREE, {});
  return {std::move(pf_paired), std::move(pf_unpaired)};
}

BoltzPseudofreeCfg::BoltzPseudofreeCfg(const PseudofreeCfg& pf)
    : paired(BuildBoltz(pf.paired)), unpaired(BuildBoltz(pf.unpaired)),
      unpaired_sum_log(BuildUnpairedSumLog(pf.unpaired)) {}

void BoltzPseudofreeCfg::Verify(const Primary& r) const {
  if (!paired.empty())
    verify(paired.size() == r.size(), "pseudofree paired must be same length as seq");
  if (!unpaired.empty())
    verify(unpaired.size() == r.size(), "pseudofree unpaired must be same length as seq");
}

}  // namespace mrna::erg
