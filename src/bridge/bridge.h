// Copyright 2016 E.
#ifndef BRIDGE_BRIDGE_H_
#define BRIDGE_BRIDGE_H_

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "common.h"
#include "model/context.h"
#include "partition/partition.h"
#include "util/argparse.h"

namespace mrna {
namespace bridge {

class RnaPackage {
 public:
  virtual ~RnaPackage() = default;

  virtual energy_t Efn(const secondary_t& secondary, std::string* desc = nullptr) const = 0;
  virtual computed_t Fold(const primary_t& r) const = 0;
  virtual int Suboptimal(
      fold::SuboptimalCallback fn, const primary_t& r, energy_t energy_delta) const = 0;
  virtual std::vector<computed_t> SuboptimalIntoVector(
      const primary_t& r, energy_t energy_delta) const = 0;
  virtual std::pair<partition::partition_t, partition::probabilities_t> Partition(
      const primary_t& r) const = 0;
};

const std::map<std::string, opt_t> BRIDGE_OPTIONS = {{"r", {"rnastructure"}}, {"k", {"memerna"}}};

std::unique_ptr<RnaPackage> RnaPackageFromArgParse(const ArgParse& args);

}  // namespace bridge
}  // namespace mrna

#endif  // BRIDGE_BRIDGE_H_
