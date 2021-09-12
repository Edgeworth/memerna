// Copyright 2016 E.
#ifndef MEMERNA_RNASTRUCTURE_H_
#define MEMERNA_RNASTRUCTURE_H_

#include "bridge/bridge.h"
#include "common.h"
#include "miles_rnastructure/include/algorithm.h"
#include "miles_rnastructure/include/alltrace.h"
#include "miles_rnastructure/include/pfunction.h"
#include "miles_rnastructure/include/rna_library.h"

namespace memerna {
namespace bridge {

class Rnastructure : public RnaPackage {
public:
  Rnastructure(const std::string& data_path, bool use_lyngso_);
  Rnastructure(const Rnastructure&) = delete;
  Rnastructure& operator=(Rnastructure&) = delete;

  virtual energy_t Efn(const secondary_t& secondary, std::string* desc = nullptr) const override;
  virtual computed_t Fold(const primary_t& r) const override;
  virtual int Suboptimal(
      fold::SuboptimalCallback fn, const primary_t& r, energy_t energy_delta) const override;
  virtual std::vector<computed_t> SuboptimalIntoVector(
      const primary_t& r, energy_t energy_delta) const override;
  virtual std::pair<partition::partition_t, partition::probabilities_t> Partition(
      const primary_t& r) const override;

  computed_t FoldAndDpTable(const primary_t& r, dp_state_t& dp_state) const;

private:
  const std::unique_ptr<datatable> data;
  const bool use_lyngso;

  std::unique_ptr<structure> LoadStructure(const primary_t& r) const;
  std::unique_ptr<structure> LoadStructure(const secondary_t& s) const;
};

}  // namespace bridge
}  // namespace memerna

#endif  // MEMERNA_RNASTRUCTURE_H_
