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

  virtual Energy Efn(const Secondary& secondary, std::string* desc = nullptr) const = 0;
  virtual Computed Fold(const Primary& r) const = 0;
  virtual int Suboptimal(
      subopt::SuboptimalCallback fn, const Primary& r, Energy energy_delta) const = 0;
  virtual std::vector<Computed> SuboptimalIntoVector(
      const Primary& r, Energy energy_delta) const = 0;
  virtual std::pair<partition::Partition, Probabilities> Partition(const Primary& r) const = 0;
};

const std::map<std::string, Opt> BRIDGE_OPTS = {{"r", {"rnastructure"}}, {"k", {"memerna"}}};

std::unique_ptr<RnaPackage> RnaPackageFromArgParse(const ArgParse& args);

}  // namespace mrna::bridge

#endif  // BRIDGE_BRIDGE_H_
