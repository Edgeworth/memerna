// Copyright 2016 Eliot Courtney.
#ifndef API_BRIDGE_MEMERNA_H_
#define API_BRIDGE_MEMERNA_H_

#include <optional>
#include <utility>
#include <vector>

#include "api/bridge/bridge.h"
#include "api/ctx/ctx.h"
#include "api/energy/energy.h"
#include "api/pfn.h"
#include "api/subopt/subopt.h"
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

  [[nodiscard]] erg::EnergyResult Efn(const Primary& r, const Secondary& s, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf, const Ctds* given_ctd = nullptr,
      bool build_structure = false) const override;

  [[nodiscard]] FoldResult Fold(const Primary& r, std::optional<MfeAlg> alg, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf, const trace::TraceCfg& trace_cfg) const override;

  [[nodiscard]] int64_t Subopt(const Primary& r, std::optional<MfeAlg> mfe_alg,
      std::optional<SuboptAlg> alg, erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf,
      const subopt::SuboptCallback& fn, subopt::SuboptCfg subopt_cfg) const override;

  [[nodiscard]] std::vector<subopt::SuboptResult> SuboptIntoVector(const Primary& r,
      std::optional<MfeAlg> mfe_alg, std::optional<SuboptAlg> alg, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf, subopt::SuboptCfg subopt_cfg) const override;

  [[nodiscard]] pfn::PfnResult Pfn(const Primary& r, std::optional<PfnAlg> alg, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf) const override;

  static Memerna FromArgParse(const ArgParse& args);

 private:
  Ctx ctx_;
};

}  // namespace mrna::bridge

#endif  // API_BRIDGE_MEMERNA_H_
