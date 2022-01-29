// Copyright 2016 E.
#ifndef BRIDGE_MEMERNA_H_
#define BRIDGE_MEMERNA_H_

#include <string>
#include <utility>
#include <vector>

#include "bridge/bridge.h"
#include "compute/energy/energy.h"
#include "compute/partition/partition.h"
#include "compute/subopt/subopt.h"
#include "context/ctx.h"
#include "model/model.h"
#include "model/secondary.h"
#include "model/primary.h"
#include "util/argparse.h"

namespace mrna::bridge {

// Note that only one energy model can be loaded at a time.
class Memerna : public RnaPackage {
 public:
  explicit Memerna(Ctx ctx) : ctx_(std::move(ctx)) {}

  Memerna(Memerna&& o) = default;
  Memerna& operator=(Memerna&&) = default;

  Memerna(const Memerna&) = delete;
  Memerna& operator=(const Memerna&) = delete;

  energy::EnergyResult Efn(Primary r, Secondary s, std::string* desc = nullptr) const override;
  FoldResult Fold(Primary r) const override;
  int Suboptimal(subopt::SuboptCallback fn, Primary r, Energy delta) const override;
  std::vector<subopt::SuboptResult> SuboptimalIntoVector(Primary r, Energy delta) const override;
  partition::PartitionResult Partition(Primary r) const override;

  static Memerna FromArgParse(const ArgParse& args);

 private:
  Ctx ctx_;
};

}  // namespace mrna::bridge

#endif  // BRIDGE_MEMERNA_H_
