// Copyright 2016 E.
#ifndef BRIDGE_RNASTRUCTURE_H_
#define BRIDGE_RNASTRUCTURE_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "bridge/bridge.h"
#include "compute/energy/energy.h"
#include "compute/partition/partition.h"
#include "compute/subopt/subopt.h"
#include "ctx/ctx.h"
#include "miles_rnastructure/include/algorithm.h"
#include "miles_rnastructure/include/alltrace.h"
#include "miles_rnastructure/include/pfunction.h"
#include "miles_rnastructure/include/rna_library.h"
#include "model/model.h"
#include "model/secondary.h"

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
  ctx::FoldResult Fold(const Primary& r) const override;
  int Suboptimal(subopt::SuboptCallback fn, const Primary& r, Energy delta) const override;
  std::vector<subopt::SuboptResult> SuboptimalIntoVector(
      const Primary& r, Energy delta) const override;
  part::PartResult Partition(const Primary& r) const override;
  // Runs the Ding & Lawrence stochastic sample algorithm. Note that the energies in SuboptResult
  // are meaningless.
  std::vector<subopt::SuboptResult> StochasticSampleIntoVector(
      const Primary& r, int num_samples) const;

  // TODO: Can be replaced by Fold now?
  ctx::FoldResult FoldAndDpTable(const Primary& r, dp_state_t* dp_state) const;

  static RNAstructure FromArgParse(const ArgParse& args);

 private:
  std::unique_ptr<datatable> data_;
  bool use_lyngso_;

  struct PartitionFunctionState {
    const PFPRECISION scaling = 1.0;  // TODO return scaling to 0.6.
    std::unique_ptr<DynProgArray<PFPRECISION>> w;
    std::unique_ptr<DynProgArray<PFPRECISION>> v;
    std::unique_ptr<DynProgArray<PFPRECISION>> wmb;
    std::unique_ptr<DynProgArray<PFPRECISION>> wl;
    std::unique_ptr<DynProgArray<PFPRECISION>> wlc;
    std::unique_ptr<DynProgArray<PFPRECISION>> wmbl;
    std::unique_ptr<DynProgArray<PFPRECISION>> wcoax;
    std::unique_ptr<PFPRECISION[]> w5;
    std::unique_ptr<PFPRECISION[]> w3;
    std::unique_ptr<pfdatatable> pfdata;
    std::unique_ptr<forceclass> fce;
    std::unique_ptr<bool[]> lfce;
    std::unique_ptr<bool[]> mod;
    PartitionFunctionState(int N, datatable* data);
  };

  PartitionFunctionState CallPartitionFunction(structure* struc) const;

  std::unique_ptr<structure> LoadStructure(const Primary& r) const;
  std::unique_ptr<structure> LoadStructure(const Primary& r, const Secondary& s) const;
};

}  // namespace mrna::bridge

#endif  // BRIDGE_RNASTRUCTURE_H_
