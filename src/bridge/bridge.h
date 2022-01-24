// Copyright 2016 Eliot Courtney.
#ifndef BRIDGE_BRIDGE_H_
#define BRIDGE_BRIDGE_H_

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "compute/boltz_dp.h"
#include "compute/partition/partition.h"
#include "compute/subopt/subopt.h"
#include "model/context.h"
#include "model/secondary.h"
#include "util/argparse.h"

namespace mrna::bridge {

class RnaPackage {
 public:
  virtual ~RnaPackage() = default;

  virtual energy::EnergyResult Efn(Primary r, Secondary s, std::string* desc = nullptr) const = 0;
  virtual FoldResult Fold(Primary r) const = 0;
  virtual int Suboptimal(subopt::SuboptCallback fn, Primary r, Energy energy_delta) const = 0;
  virtual std::vector<subopt::SuboptResult> SuboptimalIntoVector(
      Primary r, Energy energy_delta) const = 0;
  virtual partition::PartitionResult Partition(Primary r) const = 0;

  static std::unique_ptr<RnaPackage> FromArgParse(const ArgParse& args);
};

const std::map<std::string, Opt> BRIDGE_OPTS = {{"r", {"rnastructure"}}, {"k", {"memerna"}}};

}  // namespace mrna::bridge

#endif  // BRIDGE_BRIDGE_H_
