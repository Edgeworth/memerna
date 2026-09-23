// Copyright 2025 Eliot Courtney.
#ifndef BACKENDS_BASE_ENERGY_MODEL_H_
#define BACKENDS_BASE_ENERGY_MODEL_H_

#include <memory>
#include <random>
#include <string>

#include "api/ctx/backend_cfg.h"
#include "api/energy/energy.h"
#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "backends/common/base/model_base.h"
#include "backends/common/base/parse.h"
#include "backends/common/model_mixin.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"

namespace mrna::md::base {

class Model : public ModelBase, public ModelMixin<Model> {
 public:
  static constexpr auto KIND = BackendKind::BASE;

  static bool IsSupported(const BackendCfg& backend_cfg, const erg::EnergyCfg& cfg,
      const erg::PseudofreeCfg& pf, std::string* reason = nullptr);

  // ModelMixin:
  template <erg::EnergyCfg CFG, class PseudofreeCfgT>
  EnergyResult SubEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd,
      const PseudofreeCfgT& pf, int st, int en, bool build_structure = false) const;
  template <erg::EnergyCfg CFG, class PseudofreeCfgT>
  EnergyResult TotalEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd,
      const PseudofreeCfgT& pf, bool build_structure = false) const;

  bool IsValid(std::string* reason = nullptr) const { return base::ModelIsValid(*this, reason); }

  void LoadFromModelPath(const std::string& path) { base::LoadFromModelPath(*this, path); }

  void LoadRandom(const BackendCfg& /*cfg*/, const RandomModelCfg& random_cfg, std::mt19937& eng) {
    LoadRandomModel(*this, eng, random_cfg.min_energy.ToFlt<double>(),
        random_cfg.max_energy.ToFlt<double>(), RAND_MAX_HAIRPIN_SZ, RAND_MAX_NUM_HAIRPIN);
  }

 protected:
  friend class ModelMixin<Model>;

  // This is protected to prevent construction on the stack, since this structure
  // is large.
  Model() = default;

 private:
  template <erg::EnergyCfg CFG, class PseudofreeCfgT>
  Energy SubEnergyInternal(const Primary& r, const Secondary& s, const PseudofreeCfgT& pf, int st,
      int en, bool use_given_ctds, Ctds* ctd, std::unique_ptr<Structure>* struc) const;
};

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_ENERGY_MODEL_H_
