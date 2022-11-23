// Copyright 2016 Eliot Courtney.
#ifndef BRIDGE_RNASTRUCTURE_H_
#define BRIDGE_RNASTRUCTURE_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "bridge/bridge.h"
#include "compute/energy/model.h"
#include "compute/partition/partition.h"
#include "compute/subopt/subopt.h"
#include "ctx/ctx.h"
#include "model/constants.h"
#include "model/secondary.h"
#include "rnastructure_bridge/include/algorithm.h"
#include "rnastructure_bridge/include/alltrace.h"
#include "rnastructure_bridge/include/pfunction.h"
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

  energy::EnergyResult Efn(
      const Primary& r, const Secondary& s, std::string* desc = nullptr) const override;
  [[nodiscard]] ctx::FoldResult Fold(const Primary& r) const override;
  [[nodiscard]] int Suboptimal(
      subopt::SuboptCallback fn, const Primary& r, Energy delta) const override;
  [[nodiscard]] std::vector<subopt::SuboptResult> SuboptimalIntoVector(
      const Primary& r, Energy delta) const override;
  [[nodiscard]] part::PartResult Partition(const Primary& r) const override;
  // Runs the Ding & Lawrence stochastic sample algorithm. Note that the energies in SuboptResult
  // are meaningless.
  [[nodiscard]] std::vector<subopt::SuboptResult> StochasticSampleIntoVector(
      const Primary& r, int num_samples) const;

  // TODO(2): Can be replaced by Fold now?
  ctx::FoldResult FoldAndDpTable(const Primary& r, dp_state_t* dp_state) const;

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

#endif  // BRIDGE_RNASTRUCTURE_H_
