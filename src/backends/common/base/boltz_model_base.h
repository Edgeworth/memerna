// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_COMMON_BASE_BOLTZ_MODEL_BASE_H_
#define BACKENDS_COMMON_BASE_BOLTZ_MODEL_BASE_H_

#include <cassert>
#include <string>
#include <unordered_map>

#include "backends/common/base/model_base.h"
#include "model/base.h"
#include "model/energy.h"

// GCC 11 doesn't handle constructors for MPFR here well - takes a lot of time
// and memory. We set these values after initialization anyway, so skip the
// initializer for MPFR.
#ifdef MRNA_USE_MPFR
#define BOLTZ_MEMBER_INIT
#else
#define BOLTZ_MEMBER_INIT = {}
#endif

namespace mrna::md::base {

template <typename M>
class BoltzModelBase {
 public:
  BoltzModelBase() = delete;

  BoltzEnergy stack[4][4][4][4] BOLTZ_MEMBER_INIT;
  BoltzEnergy terminal[4][4][4][4] BOLTZ_MEMBER_INIT;
  BoltzEnergy internal_init[ModelBase::INITIATION_CACHE_SZ] BOLTZ_MEMBER_INIT;
  BoltzEnergy internal_1x1[4][4][4][4][4][4] BOLTZ_MEMBER_INIT;
  BoltzEnergy internal_1x2[4][4][4][4][4][4][4] BOLTZ_MEMBER_INIT;
  BoltzEnergy internal_2x2[4][4][4][4][4][4][4][4] BOLTZ_MEMBER_INIT;
  BoltzEnergy internal_2x3_mismatch[4][4][4][4] BOLTZ_MEMBER_INIT;
  BoltzEnergy internal_other_mismatch[4][4][4][4] BOLTZ_MEMBER_INIT;
  BoltzEnergy internal_asym BOLTZ_MEMBER_INIT;
  BoltzEnergy internal_au_penalty BOLTZ_MEMBER_INIT;
  BoltzEnergy internal_gu_penalty BOLTZ_MEMBER_INIT;
  BoltzEnergy bulge_init[ModelBase::INITIATION_CACHE_SZ] BOLTZ_MEMBER_INIT;
  BoltzEnergy bulge_special_c BOLTZ_MEMBER_INIT;
  BoltzEnergy hairpin_init[ModelBase::INITIATION_CACHE_SZ] BOLTZ_MEMBER_INIT;
  BoltzEnergy hairpin_uu_ga_first_mismatch BOLTZ_MEMBER_INIT;
  BoltzEnergy hairpin_gg_first_mismatch BOLTZ_MEMBER_INIT;
  BoltzEnergy hairpin_special_gu_closure BOLTZ_MEMBER_INIT;
  BoltzEnergy hairpin_c3_loop BOLTZ_MEMBER_INIT;
  BoltzEnergy hairpin_all_c_a BOLTZ_MEMBER_INIT;
  BoltzEnergy hairpin_all_c_b BOLTZ_MEMBER_INIT;
  std::unordered_map<std::string, BoltzEnergy> hairpin;
  BoltzEnergy multiloop_a BOLTZ_MEMBER_INIT;
  BoltzEnergy multiloop_b BOLTZ_MEMBER_INIT;
  BoltzEnergy multiloop_c BOLTZ_MEMBER_INIT;
  BoltzEnergy dangle5[4][4][4] BOLTZ_MEMBER_INIT;
  BoltzEnergy dangle3[4][4][4] BOLTZ_MEMBER_INIT;
  BoltzEnergy coax_mismatch_non_contiguous BOLTZ_MEMBER_INIT;
  BoltzEnergy coax_mismatch_wc_bonus BOLTZ_MEMBER_INIT;
  BoltzEnergy coax_mismatch_gu_bonus BOLTZ_MEMBER_INIT;
  BoltzEnergy au_penalty BOLTZ_MEMBER_INIT;
  BoltzEnergy gu_penalty BOLTZ_MEMBER_INIT;

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

  [[nodiscard]] const M& m() const { return m_; }

 protected:
  M m_;

  // This is private to prevent construction on the stack, since this structure
  // can be very large if arbitrary precision floats are enabled.
  explicit BoltzModelBase(const M::Ptr& m) : m_(*m) { LoadBoltzModel(*this); }
};

}  // namespace mrna::md::base

#undef BOLTZ_MEMBER_INIT

#endif  // BACKENDS_COMMON_BASE_BOLTZ_MODEL_BASE_H_
