// Copyright 2016 E.
#ifndef BRIDGE_MEMERNA_H_
#define BRIDGE_MEMERNA_H_

#include <string>
#include <utility>
#include <vector>

#include "bridge/bridge.h"
#include "model/context.h"

namespace mrna::bridge {

// Note that only one energy model can be loaded at a time.
class Memerna : public RnaPackage {
 public:
  Memerna(energy::EnergyModelPtr em_, context_opt_t options_) : em(em_), options(options_) {}

  Memerna(const Memerna&) = delete;
  Memerna& operator=(const Memerna&) = delete;

  energy_t Efn(const secondary_t& secondary, std::string* desc = nullptr) const override;
  computed_t Fold(const primary_t& r) const override;
  int Suboptimal(
      subopt::SuboptimalCallback fn, const primary_t& r, energy_t energy_delta) const override;
  std::vector<computed_t> SuboptimalIntoVector(
      const primary_t& r, energy_t energy_delta) const override;
  std::pair<partition::partition_t, partition::probabilities_t> Partition(
      const primary_t& r) const override;

 private:
  const energy::EnergyModelPtr em;
  const context_opt_t options;
};

}  // namespace mrna::bridge

#endif  // BRIDGE_MEMERNA_H_
