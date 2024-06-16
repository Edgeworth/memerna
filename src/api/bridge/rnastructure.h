// Copyright 2016 Eliot Courtney.
#ifndef API_BRIDGE_RNASTRUCTURE_H_
#define API_BRIDGE_RNASTRUCTURE_H_

#include <memory>
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

  erg::EnergyResult Efn(
      const Primary& r, const Secondary& s, std::string* desc = nullptr) const override;
  [[nodiscard]] FoldResult Fold(const Primary& r) const override;
  [[nodiscard]] int Subopt(
      subopt::SuboptCallback fn, const Primary& r, Energy delta) const override;
  [[nodiscard]] std::vector<subopt::SuboptResult> SuboptIntoVector(
      const Primary& r, Energy delta) const override;
  [[nodiscard]] pfn::PfnResult Pfn(const Primary& r) const override;
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
