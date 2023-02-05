// Copyright 2023 Eliot Courtney.
#ifndef MODELS_T04_ENERGY_BOLTZ_MIXIN_H_
#define MODELS_T04_ENERGY_BOLTZ_MIXIN_H_

#include <cassert>
#include <memory>
#include <string>
#include <unordered_map>

#include "model/base.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/structure.h"
#include "models/common/boltz.h"
#include "models/t04/energy/model.h"
#include "models/t04/energy/model_mixin.h"

namespace mrna::md::t04::erg {

template <typename T>
  requires std::is_base_of_v<T04ModelMixin, T>
class T04BoltzMixin {
 public:
  BoltzEnergy stack[4][4][4][4]{};
  BoltzEnergy terminal[4][4][4][4]{};
  BoltzEnergy internal_init[T::INITIATION_CACHE_SZ]{};
  BoltzEnergy internal_1x1[4][4][4][4][4][4]{};
  BoltzEnergy internal_1x2[4][4][4][4][4][4][4]{};
  BoltzEnergy internal_2x2[4][4][4][4][4][4][4][4]{};
  BoltzEnergy internal_2x3_mismatch[4][4][4][4]{};
  BoltzEnergy internal_other_mismatch[4][4][4][4]{};
  BoltzEnergy internal_asym{};
  BoltzEnergy internal_au_penalty{};
  BoltzEnergy internal_gu_penalty{};
  BoltzEnergy bulge_init[T::INITIATION_CACHE_SZ]{};
  BoltzEnergy bulge_special_c{};
  BoltzEnergy hairpin_init[T::INITIATION_CACHE_SZ]{};
  BoltzEnergy hairpin_uu_ga_first_mismatch{}, hairpin_gg_first_mismatch{},
      hairpin_special_gu_closure{}, hairpin_c3_loop{}, hairpin_all_c_a{}, hairpin_all_c_b{};
  std::unordered_map<std::string, BoltzEnergy> hairpin;
  BoltzEnergy multiloop_hack_a{}, multiloop_hack_b{};
  BoltzEnergy dangle5[4][4][4]{};
  BoltzEnergy dangle3[4][4][4]{};
  BoltzEnergy coax_mismatch_non_contiguous{}, coax_mismatch_wc_bonus{}, coax_mismatch_gu_bonus{};
  BoltzEnergy au_penalty{};
  BoltzEnergy gu_penalty{};

  const T& em() const { return em_; }

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
    return em_.Hairpin(r, st, en, s).Boltz();
  }

  BoltzEnergy Bulge(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const {
    return em_.Bulge(r, ost, oen, ist, ien, s).Boltz();
  }

  BoltzEnergy InternalLoop(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const {
    return em_.InternalLoop(r, ost, oen, ist, ien, s).Boltz();
  }

  BoltzEnergy TwoLoop(const Primary& r, int ost, int oen, int ist, int ien,
      std::unique_ptr<Structure>* s = nullptr) const {
    return em_.TwoLoop(r, ost, oen, ist, ien, s).Boltz();
  }

 protected:
  explicit T04BoltzMixin(const T& em_unshadow) : em_(em_unshadow) {
    // Force this to be false to not include bulge states for the partition
    // function.
    em_.cfg.bulge_states = false;

    FILL_BOLTZ(stack);
    FILL_BOLTZ(terminal);
    FILL_BOLTZ(internal_init);
    FILL_BOLTZ(internal_1x1);
    FILL_BOLTZ(internal_1x2);
    FILL_BOLTZ(internal_2x2);
    FILL_BOLTZ(internal_2x3_mismatch);
    FILL_BOLTZ(internal_other_mismatch);
    FILL_BOLTZ(internal_asym);
    FILL_BOLTZ(internal_au_penalty);
    FILL_BOLTZ(internal_gu_penalty);
    FILL_BOLTZ(bulge_init);
    FILL_BOLTZ(bulge_special_c);
    FILL_BOLTZ(hairpin_init);
    FILL_BOLTZ(hairpin_uu_ga_first_mismatch);
    FILL_BOLTZ(hairpin_gg_first_mismatch);
    FILL_BOLTZ(hairpin_special_gu_closure);
    FILL_BOLTZ(hairpin_c3_loop);
    FILL_BOLTZ(hairpin_all_c_a);
    FILL_BOLTZ(hairpin_all_c_b);
    FILL_BOLTZ(multiloop_hack_a);
    FILL_BOLTZ(multiloop_hack_b);
    FILL_BOLTZ(dangle5);
    FILL_BOLTZ(dangle3);
    FILL_BOLTZ(coax_mismatch_non_contiguous);
    FILL_BOLTZ(coax_mismatch_wc_bonus);
    FILL_BOLTZ(coax_mismatch_gu_bonus);
    FILL_BOLTZ(au_penalty);
    FILL_BOLTZ(gu_penalty);

    for (const auto& kv : em_.hairpin) hairpin[kv.first] = kv.second.Boltz();
  }

  T em_;
};

}  // namespace mrna::md::t04::erg

#endif  // MODELS_T04_ENERGY_BOLTZ_MIXIN_H_
