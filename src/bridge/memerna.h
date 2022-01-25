// Copyright 2016 E.
#ifndef BRIDGE_MEMERNA_H_
#define BRIDGE_MEMERNA_H_

#include <string>
#include <utility>
#include <vector>

#include "bridge/bridge.h"
#include "compute/energy/energy.h"
#include "compute/energy/model.h"
#include "compute/partition/partition.h"
#include "compute/subopt/subopt.h"
#include "model/config.h"
#include "model/context.h"
#include "model/model.h"
#include "model/secondary.h"

namespace mrna::bridge {

// Note that only one energy model can be loaded at a time.
class Memerna : public RnaPackage {
 public:
  Memerna(energy::EnergyModel em, ModelCfg cfg) : em_(std::move(em)), cfg_(std::move(cfg)) {}

  Memerna(const Memerna&) = delete;
  Memerna& operator=(const Memerna&) = delete;

  energy::EnergyResult Efn(Primary r, Secondary s, std::string* desc = nullptr) const override;
  FoldResult Fold(Primary r) const override;
  int Suboptimal(subopt::SuboptCallback fn, Primary r, Energy energy_delta) const override;
  std::vector<subopt::SuboptResult> SuboptimalIntoVector(
      Primary r, Energy energy_delta) const override;
  partition::PartitionResult Partition(Primary r) const override;

 private:
  energy::EnergyModel em_;
  ModelCfg cfg_;
};

}  // namespace mrna::bridge

#endif  // BRIDGE_MEMERNA_H_
