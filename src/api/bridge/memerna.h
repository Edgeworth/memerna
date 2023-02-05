// Copyright 2016 Eliot Courtney.
#ifndef BRIDGE_MEMERNA_H_
#define BRIDGE_MEMERNA_H_

#include <string>
#include <utility>
#include <vector>

#include "api/part.h"
#include "api/subopt/subopt.h"
#include "bridge/bridge.h"
#include "ctx/ctx.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/argparse.h"

namespace mrna::bridge {

// Note that only one energy model can be loaded at a time.
class Memerna : public RnaPackage {
 public:
  explicit Memerna(ctx::Ctx ctx) : ctx_(std::move(ctx)) {}
  ~Memerna() override = default;

  Memerna(Memerna&& o) = default;
  Memerna& operator=(Memerna&&) = default;

  Memerna(const Memerna&) = delete;
  Memerna& operator=(const Memerna&) = delete;

  erg::EnergyResult Efn(
      const Primary& r, const Secondary& s, std::string* desc = nullptr) const override;
  [[nodiscard]] ctx::FoldResult Fold(const Primary& r) const override;
  [[nodiscard]] int Suboptimal(
      subopt::SuboptCallback fn, const Primary& r, Energy delta) const override;
  [[nodiscard]] std::vector<subopt::SuboptResult> SuboptimalIntoVector(
      const Primary& r, Energy delta) const override;
  [[nodiscard]] part::PartResult Partition(const Primary& r) const override;

  static Memerna FromArgParse(const ArgParse& args);

 private:
  ctx::Ctx ctx_;
};

}  // namespace mrna::bridge

#endif  // BRIDGE_MEMERNA_H_
