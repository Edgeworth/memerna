// Copyright 2016 E.
#ifndef BRIDGE_RNASTRUCTURE_H_
#define BRIDGE_RNASTRUCTURE_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "bridge/bridge.h"
#include "miles_rnastructure/include/algorithm.h"
#include "miles_rnastructure/include/alltrace.h"
#include "miles_rnastructure/include/pfunction.h"
#include "miles_rnastructure/include/rna_library.h"

namespace mrna::bridge {

class RNAstructure : public RnaPackage {
 public:
  RNAstructure(const std::string& data_path, bool use_lyngso);
  RNAstructure(RNAstructure&&) = default;
  RNAstructure& operator=(RNAstructure&&) = default;

  RNAstructure(const RNAstructure&) = delete;
  RNAstructure& operator=(const RNAstructure&) = delete;

  energy::EnergyResult Efn(
      const Primary& r, const Secondary& s, std::string* desc = nullptr) const override;
  FoldResult Fold(const Primary& r) const override;
  int Suboptimal(subopt::SuboptCallback fn, const Primary& r, Energy energy_delta) const override;
  std::vector<subopt::SuboptResult> SuboptimalIntoVector(
      const Primary& r, Energy energy_delta) const override;
  partition::PartitionResult Partition(const Primary& r) const override;

  // TODO: Can be replaced by Fold now?
  FoldResult FoldAndDpTable(const Primary& r, dp_state_t* dp_state) const;

 private:
  std::unique_ptr<datatable> data_;
  bool use_lyngso_;

  std::unique_ptr<structure> LoadStructure(const Primary& r) const;
  std::unique_ptr<structure> LoadStructure(const Primary& r, const Secondary& s) const;
};

}  // namespace mrna::bridge

#endif  // BRIDGE_RNASTRUCTURE_H_
