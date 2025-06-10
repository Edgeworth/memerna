// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_COMMON_BASE_BOLTZ_MODEL_BASE_H_
#define BACKENDS_COMMON_BASE_BOLTZ_MODEL_BASE_H_

#include <cassert>
#include <memory>
#include <string>
#include <unordered_map>

#include "backends/common/base/model_base.h"
#include "model/base.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/structure.h"

namespace mrna::md::base {

template <typename M>
class BoltzModelBase {
 public:
  BoltzModelBase() = delete;

  BoltzEnergy stack[4][4][4][4]{};
  BoltzEnergy terminal[4][4][4][4]{};
  BoltzEnergy internal_init[ModelBase::INITIATION_CACHE_SZ]{};
  BoltzEnergy internal_1x1[4][4][4][4][4][4]{};
  BoltzEnergy internal_1x2[4][4][4][4][4][4][4]{};
  BoltzEnergy internal_2x2[4][4][4][4][4][4][4][4]{};
  BoltzEnergy internal_2x3_mismatch[4][4][4][4]{};
  BoltzEnergy internal_other_mismatch[4][4][4][4]{};
  BoltzEnergy internal_asym{};
  BoltzEnergy internal_au_penalty{};
  BoltzEnergy internal_gu_penalty{};
  BoltzEnergy bulge_init[ModelBase::INITIATION_CACHE_SZ]{};
  BoltzEnergy bulge_special_c{};
  BoltzEnergy hairpin_init[ModelBase::INITIATION_CACHE_SZ]{};
  BoltzEnergy hairpin_uu_ga_first_mismatch{};
  BoltzEnergy hairpin_gg_first_mismatch{};
  BoltzEnergy hairpin_special_gu_closure{};
  BoltzEnergy hairpin_c3_loop{};
  BoltzEnergy hairpin_all_c_a{};
  BoltzEnergy hairpin_all_c_b{};
  std::unordered_map<std::string, BoltzEnergy> hairpin;
  BoltzEnergy multiloop_a{};
  BoltzEnergy multiloop_b{};
  BoltzEnergy multiloop_c{};
  BoltzEnergy dangle5[4][4][4]{};
  BoltzEnergy dangle3[4][4][4]{};
  BoltzEnergy coax_mismatch_non_contiguous{};
  BoltzEnergy coax_mismatch_wc_bonus{};
  BoltzEnergy coax_mismatch_gu_bonus{};
  BoltzEnergy au_penalty{};
  BoltzEnergy gu_penalty{};

  BoltzEnergy InternalLoopAuGuPenalty(Base stb, Base enb) const {
    assert(IsBase(stb) && IsBase(enb));
    if (IsAuPair(stb, enb)) return internal_au_penalty;
    if (IsGuPair(stb, enb)) return internal_gu_penalty;
    return 1.0;
  }

  BoltzEnergy AuGuPenalty(Base stb, Base enb) const {
    assert(IsBase(stb) && IsBase(enb));
    if (IsAuPair(stb, enb)) return au_penalty;
    if (IsGuPair(stb, enb)) return gu_penalty;
    return 1.0;
  }

  BoltzEnergy MismatchCoaxial(
      Base five_top, Base mismatch_top, Base mismatch_bot, Base three_bot) const {
    assert(IsBase(five_top) && IsBase(mismatch_top) && IsBase(mismatch_bot) && IsBase(three_bot));
    BoltzEnergy coax =
        terminal[five_top][mismatch_top][mismatch_bot][three_bot] * coax_mismatch_non_contiguous;
    if (IsWcPair(mismatch_top, mismatch_bot))
      coax *= coax_mismatch_wc_bonus;
    else if (IsGuPair(mismatch_top, mismatch_bot))
      coax *= coax_mismatch_gu_bonus;
    return coax;
  }

  BoltzEnergy Hairpin(
      const Primary& r, int st, int en, std::unique_ptr<Structure>* s = nullptr) const {
    return m_.Hairpin(r, st, en, s).Boltz();
  }

  BoltzEnergy Bulge(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const {
    return m_.Bulge(r, ost, oen, ist, ien, s).Boltz();
  }

  BoltzEnergy InternalLoop(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const {
    return m_.InternalLoop(r, ost, oen, ist, ien, s).Boltz();
  }

  BoltzEnergy TwoLoop(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const {
    return m_.TwoLoop(r, ost, oen, ist, ien, s).Boltz();
  }

  [[nodiscard]] const M& m() const { return m_; }

 protected:
  M m_;

  // This is private to prevent construction on the stack, since this structure
  // can be very large if arbitrary precision floats are enabled.
  explicit BoltzModelBase(const M::Ptr& m) : m_(*m) {
    // Force this to be false to not include bulge states for the partition
    // function.
    auto cfg = m_.cfg();
    cfg.bulge_states = false;
    m_.SetEnergyCfg(cfg);
    LoadBoltzModel(*this);
  }
};

}  // namespace mrna::md::base

#endif  // BACKENDS_COMMON_BASE_BOLTZ_MODEL_BASE_H_
