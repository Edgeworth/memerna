// Copyright 2016 Eliot Courtney.
#ifndef BRIDGE_MEMERNA_H_
#define BRIDGE_MEMERNA_H_

#include <string>
#include <utility>
#include <vector>

#include "bridge/bridge.h"
#include "model/config.h"
#include "model/context.h"

namespace mrna::bridge {

// Note that only one energy model can be loaded at a time.
class Memerna : public RnaPackage {
 public:
  Memerna(energy::EnergyModel em, ModelCfg cfg) : em_(std::move(em)), cfg_(std::move(cfg)) {}

  Memerna(const Memerna&) = delete;
  Memerna& operator=(const Memerna&) = delete;

  Energy Efn(const Secondary& secondary, std::string* desc = nullptr) const override;
  Computed Fold(const Primary& r) const override;
  int Suboptimal(
      subopt::SuboptimalCallback fn, const Primary& r, Energy energy_delta) const override;
  std::vector<Computed> SuboptimalIntoVector(const Primary& r, Energy energy_delta) const override;
  std::pair<partition::Partition, partition::Probabilities> Partition(
      const Primary& r) const override;

 private:
  energy::EnergyModel em_;
  ModelCfg cfg_;
};

}  // namespace mrna::bridge

#endif  // BRIDGE_MEMERNA_H_
