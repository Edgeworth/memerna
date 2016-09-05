#include "fold/fold_internal.h"
#include "parsing.h"

namespace memerna {
namespace fold {
namespace internal {

using namespace constants;
using namespace energy;

namespace {

energy_t MinEnergy(const energy_t* energy, std::size_t size) {
  energy_t min = energy[0];
  for (int i = 0; i < int(size / sizeof(energy_t)); ++i)
    min = std::min(min, energy[i]);
  return min;
}

}

int MaxNumContiguous(const primary_t& r) {
  energy_t num_contig = 0;
  energy_t max_num_contig = 0;
  base_t prev = -1;
  for (auto b : r) {
    if (b == prev) num_contig++;
    else num_contig = 1;
    prev = b;
    max_num_contig = std::max(max_num_contig, num_contig);
  }
  return max_num_contig;
}

precomp_t PrecomputeData(const primary_t& r, const energy::EnergyModel& em) {
  assert(r.size() > 0);
  precomp_t pc;
  // Initialise fast AUGU branch table
  for (base_t i = 0; i < 4; ++i)
    for (base_t j = 0; j < 4; ++j)
      pc.augubranch[i][j] = em.multiloop_hack_b + em.AuGuPenalty(i, j);

  const auto min_stack = MinEnergy(&em.stack[0][0][0][0], sizeof(em.stack));

  // Non continuous (-2.1), -4 for WC, -16 for terminal mismatch.
  pc.min_mismatch_coax = em.coax_mismatch_non_contiguous +
      std::min(std::min(em.coax_mismatch_gu_bonus, em.coax_mismatch_wc_bonus), 0) +
      MinEnergy(&em.terminal[0][0][0][0], sizeof(em.terminal));
  // Minimum of all stacking params.
  pc.min_flush_coax = min_stack;

  energy_t min_internal = MinEnergy(&em.internal_1x1[0][0][0][0][0][0], sizeof(em.internal_1x1));
  min_internal = std::min(min_internal,
      MinEnergy(&em.internal_1x2[0][0][0][0][0][0][0], sizeof(em.internal_1x2)));
  min_internal = std::min(min_internal,
      MinEnergy(&em.internal_2x2[0][0][0][0][0][0][0][0], sizeof(em.internal_2x2)));
  verify_expr(em.internal_asym >= 0,
      "min_internal optimisation does not work for negative asymmetry penalties");
  const auto min_mismatch = 2 * std::min(
      MinEnergy(&em.internal_2x3_mismatch[0][0][0][0], sizeof(em.internal_2x3_mismatch)),
      MinEnergy(&em.internal_other_mismatch[0][0][0][0], sizeof(em.internal_other_mismatch)));
  const auto min_internal_init = MinEnergy(
      &em.internal_init[4], sizeof(em.internal_init) - 4 * sizeof(em.internal_init[0]));
  min_internal = std::min(min_internal,
      min_internal_init + std::min(2 * em.internal_augu_penalty, 0) + min_mismatch);

  const auto min_bulge_init = MinEnergy(&em.bulge_init[1], sizeof(em.bulge_init) - sizeof(em.bulge_init[0]));

  energy_t states_bonus = -energy_t(round(10.0 * R * T * log(MaxNumContiguous(r))));
  energy_t min_bulge = min_bulge_init + std::min(2 * em.augu_penalty, 0) +
      min_stack + std::min(em.bulge_special_c, 0) + states_bonus;
  pc.min_twoloop_not_stack = std::min(min_bulge, min_internal);

  pc.hairpin.resize(r.size());
  std::string rna_str = parsing::PrimaryToString(r);
  for (const auto& hairpinpair : em.hairpin) {
    const auto& str = hairpinpair.first;
    verify_expr(str.size() - 2 <= hairpin_precomp_t::MAX_SPECIAL_HAIRPIN_SZ, "need to increase MAX_SPECIAL_HAIRPIN_SZ");
    auto pos = rna_str.find(str, 0);
    while (pos != std::string::npos) {
      pc.hairpin[pos].special[str.size() - 2] = hairpinpair.second;
      pos = rna_str.find(str, pos + 1);
    }
  }
  const int N = int(r.size());
  pc.hairpin[N - 1].num_c = int(r[N - 1] == C);
  for (int i = N - 2; i >= 0; --i) {
    if (r[i] == C)
      pc.hairpin[i].num_c = pc.hairpin[i + 1].num_c + 1;
  }

  return pc;
}

}
}
}
