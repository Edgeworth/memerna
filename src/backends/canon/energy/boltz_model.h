// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_CANON_ENERGY_BOLTZ_MODEL_H_
#define BACKENDS_CANON_ENERGY_BOLTZ_MODEL_H_

#include <memory>

#include "api/energy/energy_cfg.h"
#include "backends/canon/energy/model.h"
#include "backends/common/base/boltz_model_base.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/structure.h"

namespace mrna::md::canon {

class BoltzModel : public base::BoltzModelBase<Model> {
 public:
  using Ptr = std::shared_ptr<BoltzModel>;

  BoltzModel() = delete;

  static Ptr Create(const Model::Ptr& m) { return Ptr(new BoltzModel(m)); }

  template <erg::EnergyCfg CFG, class PseudofreeCfgT>
  BoltzEnergy Hairpin(const Primary& r, const PseudofreeCfgT& pf, int st, int en,
      std::unique_ptr<Structure>* s = nullptr) const {
    return Boltz(m_.Hairpin<CFG>(r, pf, st, en, s));
  }

  template <erg::EnergyCfg CFG, class PseudofreeCfgT>
  BoltzEnergy Bulge(const Primary& r, const PseudofreeCfgT& pf, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const {
    return Boltz(m_.Bulge<CFG>(r, pf, ost, oen, ist, ien, s));
  }

  template <erg::EnergyCfg CFG, class PseudofreeCfgT>
  BoltzEnergy InternalLoop(const Primary& r, const PseudofreeCfgT& pf, int ost, int oen, int ist,
      int ien, std::unique_ptr<Structure>* s = nullptr) const {
    return Boltz(m_.InternalLoop<CFG>(r, pf, ost, oen, ist, ien, s));
  }

  template <erg::EnergyCfg CFG, class PseudofreeCfgT>
  BoltzEnergy TwoLoop(const Primary& r, const PseudofreeCfgT& pf, int ost, int oen, int ist,
      int ien, std::unique_ptr<Structure>* s = nullptr) const {
    return Boltz(m_.TwoLoop<CFG>(r, pf, ost, oen, ist, ien, s));
  }

 private:
  // This is private to prevent construction on the stack, since this structure is
  // large.
  explicit BoltzModel(const Model::Ptr& m) : BoltzModelBase(m) {}
};

}  // namespace mrna::md::canon

#endif  // BACKENDS_CANON_ENERGY_BOLTZ_MODEL_H_
