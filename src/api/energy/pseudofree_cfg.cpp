// Copyright 2025 Eliot Courtney.
#include "api/energy/pseudofree_cfg.h"

#include <utility>
#include <vector>

#include "util/error.h"

namespace mrna::erg {

void RegisterOptsPseudofree(ArgParse* args) {
  args->RegisterOpt(OPT_PAIRED_PSEUDOFREE);
  args->RegisterOpt(OPT_UNPAIRED_PSEUDOFREE);
}

void PseudofreeCfg::Load(std::vector<Energy> pf_paired, std::vector<Energy> pf_unpaired) {
  paired = std::move(pf_paired);
  unpaired = std::move(pf_unpaired);
  if (!unpaired.empty()) {
    unpaired_cum.resize(unpaired.size() + 1);
    unpaired_cum[0] = ZERO_E;
    for (size_t i = 0; i < unpaired.size(); ++i)
      unpaired_cum[i + 1] = unpaired_cum[i] + unpaired[i];
  }
}

void PseudofreeCfg::Verify(const Primary& r) const {
  if (!paired.empty())
    verify(paired.size() == r.size(), "pseudofree paired must be same length as seq");
  if (!unpaired.empty())
    verify(unpaired.size() == r.size(), "pseudofree unpaired must be same length as seq");
}

PseudofreeCfg PseudofreeCfg::FromArgParse(const ArgParse& args) {
  PseudofreeCfg cfg;
  auto pf_paired = args.GetMultipleOr<Energy>(OPT_PAIRED_PSEUDOFREE, {});
  auto pf_unpaired = args.GetMultipleOr<Energy>(OPT_UNPAIRED_PSEUDOFREE, {});
  cfg.Load(std::move(pf_paired), std::move(pf_unpaired));
  return cfg;
}

void BoltzPseudofreeCfg::Load(std::vector<Energy> pf_paired, std::vector<Energy> pf_unpaired) {
  paired.resize(pf_paired.size());
  for (size_t i = 0; i < pf_paired.size(); ++i) paired[i] = pf_paired[i].Boltz();

  unpaired.resize(pf_unpaired.size());
  for (size_t i = 0; i < pf_unpaired.size(); ++i) unpaired[i] = pf_unpaired[i].Boltz();

  if (!unpaired.empty()) {
    unpaired_cum_log.resize(unpaired.size() + 1);
    unpaired_cum_log[0] = ZERO_B;
    for (size_t i = 0; i < unpaired.size(); ++i)
      unpaired_cum_log[i + 1] = unpaired_cum_log[i] + pf_unpaired[i].LogBoltz();
  }
}

void BoltzPseudofreeCfg::Verify(const Primary& r) const {
  if (!paired.empty())
    verify(paired.size() == r.size(), "pseudofree paired must be same length as seq");
  if (!unpaired.empty())
    verify(unpaired.size() == r.size(), "pseudofree unpaired must be same length as seq");
}

}  // namespace mrna::erg
