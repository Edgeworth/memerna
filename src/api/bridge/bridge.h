// Copyright 2016 Eliot Courtney.
#ifndef API_BRIDGE_BRIDGE_H_
#define API_BRIDGE_BRIDGE_H_

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "api/ctx/algorithm.h"
#include "api/ctx/ctx.h"
#include "api/energy/energy.h"
#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "api/pfn.h"
#include "api/subopt/subopt.h"
#include "api/subopt/subopt_cfg.h"
#include "api/trace/trace_cfg.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/argparse.h"

namespace mrna::bridge {

inline const Opt OPT_USE_RNASTRUCTURE =
    Opt(Opt::FLAG).LongName("rnastructure").ShortName("r").Default(false);
inline const Opt OPT_RNASTRUCTURE_DATA =
    Opt(Opt::ARG).LongName("rnastructure-data").ShortName("rd").Help("data path for RNAstructure");
inline const Opt OPT_USE_MEMERNA = Opt(Opt::FLAG).LongName("memerna").ShortName("m").Default(false);

void RegisterOpts(ArgParse* args);

class RnaPackage {
 public:
  RnaPackage() = default;
  virtual ~RnaPackage() = default;

  RnaPackage(RnaPackage&& o) = default;
  RnaPackage& operator=(RnaPackage&&) = default;

  RnaPackage(const RnaPackage&) = delete;
  RnaPackage& operator=(const RnaPackage&) = delete;

  virtual erg::EnergyResult Efn(const Primary& r, const Secondary& s, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf, const Ctds* given_ctd = nullptr,
      bool build_structure = false) const = 0;

  [[nodiscard]] virtual FoldResult Fold(const Primary& r, std::optional<MfeAlg> alg,
      erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf, const trace::TraceCfg& trace_cfg) const = 0;

  [[nodiscard]] virtual int Subopt(const Primary& r, std::optional<MfeAlg> mfe_alg,
      std::optional<SuboptAlg> alg, erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf,
      const subopt::SuboptCallback& fn, subopt::SuboptCfg subopt_cfg) const = 0;

  [[nodiscard]] virtual std::vector<subopt::SuboptResult> SuboptIntoVector(const Primary& r,
      std::optional<MfeAlg> mfe_alg, std::optional<SuboptAlg> alg, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf, subopt::SuboptCfg subopt_cfg) const = 0;

  [[nodiscard]] virtual pfn::PfnResult Pfn(const Primary& r, std::optional<PfnAlg> alg,
      erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf) const = 0;

  static std::unique_ptr<RnaPackage> FromArgParse(const ArgParse& args);
};

}  // namespace mrna::bridge

#endif  // API_BRIDGE_BRIDGE_H_
