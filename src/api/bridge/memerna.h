// Copyright 2016 E.
#ifndef API_BRIDGE_MEMERNA_H_
#define API_BRIDGE_MEMERNA_H_

#include <string>
#include <utility>
#include <vector>

#include "api/bridge/bridge.h"
#include "api/ctx/ctx.h"
#include "api/energy/energy.h"
#include "api/pfn.h"
#include "api/subopt/subopt.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/argparse.h"

namespace mrna::bridge {

// Note that only one energy model can be loaded at a time.
class Memerna : public RnaPackage {
 public:
  explicit Memerna(Ctx ctx) : ctx_(std::move(ctx)) {}
  ~Memerna() override = default;

  Memerna(Memerna&& o) = default;
  Memerna& operator=(Memerna&&) = default;

  Memerna(const Memerna&) = delete;
  Memerna& operator=(const Memerna&) = delete;

  erg::EnergyResult Efn(
      const Primary& r, const Secondary& s, std::string* desc = nullptr) const override;
  [[nodiscard]] FoldResult Fold(const Primary& r) const override;
  [[nodiscard]] int Subopt(
      subopt::SuboptCallback fn, const Primary& r, Energy delta) const override;
  [[nodiscard]] std::vector<subopt::SuboptResult> SuboptIntoVector(
      const Primary& r, Energy delta) const override;
  [[nodiscard]] pfn::PfnResult Pfn(const Primary& r) const override;

  static Memerna FromArgParse(const ArgParse& args);

 private:
  Ctx ctx_;
};

}  // namespace mrna::bridge

#endif  // API_BRIDGE_MEMERNA_H_
