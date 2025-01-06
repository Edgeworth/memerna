// Copyright 2022 Eliot Courtney.
#ifndef BACKENDS_COMMON_BASE_PRECOMP_BASE_H_
#define BACKENDS_COMMON_BASE_PRECOMP_BASE_H_
#include <algorithm>
#include <string>
#include <vector>

#include "backends/common/energy.h"
#include "model/base.h"
#include "model/branch.h"
#include "model/energy.h"
#include "model/primary.h"
#include "util/error.h"

namespace mrna::md::base {

inline constexpr int MAX_SPECIAL_HAIRPIN_SZ = 6;

// This is templated because the partition function wants to use it with a different type.
template <typename T>
struct HairpinPrecomp {
  explicit HairpinPrecomp(T init) {
    std::fill(special, special + sizeof(special) / sizeof(special[0]), init);
  }

  T special[MAX_SPECIAL_HAIRPIN_SZ + 1];
  int num_c{0};
};

template <typename HP, bool is_boltz, typename M>
std::vector<HP> PrecomputeHairpin(const Primary& r, const M& m, auto init) {
  std::vector<HP> pc(r.size(), HP(init));
  std::string rna_str = r.ToSeq();
  for (const auto& hairpinpair : m.hairpin) {
    const auto& str = hairpinpair.first;
    verify(str.size() - 2 <= MAX_SPECIAL_HAIRPIN_SZ, "need to increase MAX_SPECIAL_HAIRPIN_SZ");
    auto energy = hairpinpair.second;
    auto augu = m.AuGuPenalty(CharToBase(*str.begin()).value(), CharToBase(*str.rbegin()).value());
    if constexpr (is_boltz)
      energy *= augu;
    else
      energy += augu;

    auto pos = rna_str.find(str, 0);
    while (pos != std::string::npos) {
      pc[pos].special[str.size() - 2] = energy;
      pos = rna_str.find(str, pos + 1);
    }
  }
  const int N = static_cast<int>(r.size());
  pc[N - 1].num_c = static_cast<int>(r[N - 1] == C);
  for (int i = N - 2; i >= 0; --i)
    if (r[i] == C) pc[i].num_c = pc[i + 1].num_c + 1;
  return pc;
}

template <typename M>
class PrecompBase {
 public:
  Energy augubranch[4][4]{};
  Energy min_mismatch_coax{};
  Energy min_flush_coax{};
  Energy min_twoloop_not_stack{};

  std::vector<HairpinPrecomp<Energy>> hairpin;

  PrecompBase(Primary r, M::Ptr m) : r_(std::move(r)), m_(std::move(m)) { PrecomputeData(); }

 protected:
  Primary r_;
  M::Ptr m_;

 private:
  void PrecomputeData() {
    // Initialise fast AUGU branch table
    for (Base i = 0; i < 4; ++i)
      for (Base j = 0; j < 4; ++j) augubranch[i][j] = m_->multiloop_hack_b + m_->AuGuPenalty(i, j);

    const auto min_stack = MinEnergy(&m_->stack[0][0][0][0], sizeof(m_->stack));

    // Non continuous (-2.1), -4 for WC, -16 for terminal mismatch.
    min_mismatch_coax = m_->coax_mismatch_non_contiguous +
        std::min(std::min(m_->coax_mismatch_gu_bonus, m_->coax_mismatch_wc_bonus), ZERO_E) +
        MinEnergy(&m_->terminal[0][0][0][0], sizeof(m_->terminal));
    // Minimum of all stacking params.
    min_flush_coax = min_stack;

    Energy min_internal = MinEnergy(&m_->internal_1x1[0][0][0][0][0][0], sizeof(m_->internal_1x1));
    min_internal = std::min(
        min_internal, MinEnergy(&m_->internal_1x2[0][0][0][0][0][0][0], sizeof(m_->internal_1x2)));
    min_internal = std::min(min_internal,
        MinEnergy(&m_->internal_2x2[0][0][0][0][0][0][0][0], sizeof(m_->internal_2x2)));

    verify(m_->internal_asym >= ZERO_E,
        "min_internal optimisation does not work for negative asymmetry penalties");
    const auto min_mismatch = 2 *
        std::min(
            MinEnergy(&m_->internal_2x3_mismatch[0][0][0][0], sizeof(m_->internal_2x3_mismatch)),
            MinEnergy(
                &m_->internal_other_mismatch[0][0][0][0], sizeof(m_->internal_other_mismatch)));
    const auto min_internal_init = MinEnergy(
        &m_->internal_init[4], sizeof(m_->internal_init) - 4 * sizeof(m_->internal_init[0]));

    // Current implementation needs regular AU/GU for special internal loops.
    const auto min_internal_penalty =
        std::min(2 * std::min(m_->internal_au_penalty, m_->internal_gu_penalty), ZERO_E);
    const auto min_penalty = std::min(2 * std::min(m_->au_penalty, m_->gu_penalty), ZERO_E);

    min_internal = std::min(min_internal + min_penalty,
        min_internal_init + min_penalty + min_internal_penalty + min_mismatch);

    const auto min_bulge_init =
        MinEnergy(&m_->bulge_init[1], sizeof(m_->bulge_init) - sizeof(m_->bulge_init[0]));

    const Energy states_bonus = -E(R * T * log(MaxNumContiguous(r_)));
    const Energy min_bulge = min_bulge_init + min_penalty + min_stack +
        std::min(m_->bulge_special_c, ZERO_E) + states_bonus;
    min_twoloop_not_stack = std::min(min_bulge, min_internal);

    hairpin = PrecomputeHairpin<HairpinPrecomp<Energy>, /*is_boltz=*/false>(r_, *m_, MAX_E);
  }
};

}  // namespace mrna::md::base

#endif  // BACKENDS_COMMON_BASE_PRECOMP_BASE_H_
