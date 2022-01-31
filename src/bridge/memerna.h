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
#include "ctx/ctx.h"
#include "model/model.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/argparse.h"

namespace mrna::bridge {

// Note that only one energy model can be loaded at a time.
class Memerna : public RnaPackage {
 public:
  explicit Memerna(ctx::Ctx ctx) : ctx_(std::move(ctx)) {}

  Memerna(Memerna&& o) = default;
  Memerna& operator=(Memerna&&) = default;

  Memerna(const Memerna&) = delete;
  Memerna& operator=(const Memerna&) = delete;

  energy::EnergyResult Efn(
      const Primary& r, const Secondary& s, std::string* desc = nullptr) const override;
  ctx::FoldResult Fold(const Primary& r) const override;
  int Suboptimal(subopt::SuboptCallback fn, const Primary& r, Energy delta) const override;
  std::vector<subopt::SuboptResult> SuboptimalIntoVector(
      const Primary& r, Energy delta) const override;
  part::PartResult Partition(const Primary& r) const override;

  static Memerna FromArgParse(const ArgParse& args);

 private:
  ctx::Ctx ctx_;
};

}  // namespace mrna::bridge

#endif  // BRIDGE_MEMERNA_H_
