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
  RNAstructure(const std::string& data_path, bool use_lyngso_);
  RNAstructure(RNAstructure&&) = default;
  RNAstructure& operator=(RNAstructure&&) = default;

  RNAstructure(const RNAstructure&) = delete;
  RNAstructure& operator=(const RNAstructure&) = delete;

  energy_t Efn(const secondary_t& secondary, std::string* desc = nullptr) const override;
  computed_t Fold(const primary_t& r) const override;
  int Suboptimal(
      subopt::SuboptimalCallback fn, const primary_t& r, energy_t energy_delta) const override;
  std::vector<computed_t> SuboptimalIntoVector(
      const primary_t& r, energy_t energy_delta) const override;
  std::pair<partition::partition_t, partition::probabilities_t> Partition(
      const primary_t& r) const override;

  computed_t FoldAndDpTable(const primary_t& r, dp_state_t* dp_state) const;

 private:
  std::unique_ptr<datatable> data;
  bool use_lyngso;

  std::unique_ptr<structure> LoadStructure(const primary_t& r) const;
  std::unique_ptr<structure> LoadStructure(const secondary_t& s) const;
};

}  // namespace mrna::bridge

#endif  // BRIDGE_RNASTRUCTURE_H_
