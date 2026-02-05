// Copyright 2025 Eliot Courtney.
#ifndef BACKENDS_BASE_ENERGY_MODEL_H_
#define BACKENDS_BASE_ENERGY_MODEL_H_

#include <deque>
#include <memory>
#include <string>

#include "api/ctx/backend_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "backends/common/base/model_base.h"
#include "backends/common/base/parse.h"
#include "backends/common/model_mixin.h"

namespace mrna::md::base {

class Model : public ModelBase, public ModelMixin<Model> {
 public:
  static constexpr auto KIND = BackendKind::BASE;

  static bool IsSupported(const BackendCfg& backend_cfg, const erg::EnergyCfg& cfg,
      const erg::PseudofreeCfg& pf, std::string* reason = nullptr);

  Energy Hairpin(const Primary& r, const erg::PseudofreeCfg& pf, int st, int en,
      std::unique_ptr<Structure>* s = nullptr) const;
  Energy Bulge(const Primary& r, erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf, int ost, int oen,
      int ist, int ien, std::unique_ptr<Structure>* s = nullptr) const;
  Energy InternalLoop(const Primary& r, const erg::PseudofreeCfg& pf, int ost, int oen, int ist,
      int ien, std::unique_ptr<Structure>* s = nullptr) const;
  Energy TwoLoop(const Primary& r, erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf, int ost,
      int oen, int ist, int ien, std::unique_ptr<Structure>* s = nullptr) const;
  Energy MultiloopEnergy(const Primary& r, const Secondary& s, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf, int st, int en, std::deque<int>* branches, bool use_given_ctds,
      Ctds* ctd, std::unique_ptr<Structure>* sstruc = nullptr) const;

  // ModelMixin:
  EnergyResult SubEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd,
      erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf, int st, int en,
      bool build_structure = false) const;
  EnergyResult TotalEnergy(const Primary& r, const Secondary& s, const Ctds* given_ctd,
      erg::EnergyCfg cfg, const erg::PseudofreeCfg& pf, bool build_structure = false) const;

  bool IsValid(std::string* reason = nullptr) const { return base::ModelIsValid(*this, reason); }

  void LoadFromModelPath(const std::string& path) { base::LoadFromModelPath(*this, path); }

  void LoadRandom(std::mt19937& eng) {
    LoadRandomModel(
        *this, eng, RAND_MIN_ENERGY, RAND_MAX_ENERGY, RAND_MAX_HAIRPIN_SZ, RAND_MAX_NUM_HAIRPIN);
  }

 private:
  friend class ModelMixin<Model>;

  // This is private to prevent construction on the stack, since this structure is large.
  Model() = default;

  Energy SubEnergyInternal(const Primary& r, const Secondary& s, erg::EnergyCfg cfg,
      const erg::PseudofreeCfg& pf, int st, int en, bool use_given_ctds, Ctds* ctd,
      std::unique_ptr<Structure>* struc) const;
};

}  // namespace mrna::md::base

#endif  // BACKENDS_BASE_ENERGY_MODEL_H_
