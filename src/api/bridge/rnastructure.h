// Copyright 2016 Eliot Courtney.
#ifndef API_BRIDGE_RNASTRUCTURE_H_
#define API_BRIDGE_RNASTRUCTURE_H_

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "api/bridge/bridge.h"
#include "api/ctx/ctx.h"
#include "api/pfn.h"
#include "api/subopt/subopt.h"
#include "model/secondary.h"
#include "rnastructure_bridge/include/algorithm.h"
#include "rnastructure_bridge/include/rna_library.h"

namespace mrna::bridge {

class RNAstructure : public RnaPackage {
 public:
  RNAstructure(const std::string& data_path, bool use_lyngso);
  ~RNAstructure() override = default;

  RNAstructure(RNAstructure&&) = default;
  RNAstructure& operator=(RNAstructure&&) = default;

  RNAstructure(const RNAstructure&) = delete;
  RNAstructure& operator=(const RNAstructure&) = delete;

  erg::EnergyResult Efn(const Primary& r, const Secondary& s, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf, const Ctds* given_ctd = nullptr,
      bool build_structure = false) const override;

  [[nodiscard]] FoldResult Fold(const Primary& r, std::optional<MfeAlg> alg, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf, const trace::TraceCfg& trace_cfg) const override;

  [[nodiscard]] int Subopt(const Primary& r, std::optional<MfeAlg> mfe_alg,
      std::optional<SuboptAlg> alg, erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf,
      const subopt::SuboptCallback& fn, subopt::SuboptCfg subopt_cfg) const override;

  [[nodiscard]] std::vector<subopt::SuboptResult> SuboptIntoVector(const Primary& r,
      std::optional<MfeAlg> mfe_alg, std::optional<SuboptAlg> alg, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf, subopt::SuboptCfg subopt_cfg) const override;

  [[nodiscard]] pfn::PfnResult Pfn(const Primary& r, std::optional<PfnAlg> alg, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf) const override;

  // Runs the Ding & Lawrence stochastic sample algorithm. Note that the energies in SuboptResult
  // are meaningless.
  [[nodiscard]] std::vector<subopt::SuboptResult> StochasticSampleIntoVector(
      const Primary& r, int num_samples) const;

  // TODO(2): Can be replaced by Fold now?
  FoldResult FoldAndDpTable(const Primary& r, dp_state_t* dp_state) const;

  static RNAstructure FromArgParse(const ArgParse& args);

  static Energy ToEnergy(int energy) { return Energy::FromRaw(energy * Energy::FACTOR / 10); }

  static int16_t FromEnergy(Energy energy) {
    auto rstr_energy = energy.v * 10 / Energy::FACTOR;
    verify(int16_t(rstr_energy) == rstr_energy, "energy too big");
    return int16_t(rstr_energy);
  }

 private:
  std::unique_ptr<datatable> data_{};
  bool use_lyngso_;

  [[nodiscard]] std::unique_ptr<structure> LoadStructure(const Primary& r) const;
  [[nodiscard]] std::unique_ptr<structure> LoadStructure(
      const Primary& r, const Secondary& s) const;
};

}  // namespace mrna::bridge

#endif  // API_BRIDGE_RNASTRUCTURE_H_
