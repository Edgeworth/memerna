#include <stack>
#include "parsing.h"
#include "fold/fold.h"
#include "fold/traceback.h"
#include "fold/fold_globals.h"

namespace memerna {
namespace fold {

using constants::MAX_E;

template<typename T>
folded_rna_t FoldInternal(const rna_t& rna, fold_state_t* fold_state, T ComputeTables) {
  SetRna(rna);
  InitFold();
  auto arr = ComputeTables();
  traceback_stack_t q;
  auto energy = TraceExterior(arr, q);
  TraceStructure(arr, q);
  if (fold_state)
    fold_state->dp_table = std::move(arr);
  return {rna, p, energy};
}

folded_rna_t Fold0(const rna_t& rna, fold_state_t* fold_state) {
  return FoldInternal(rna, fold_state, ComputeTables0);
}

folded_rna_t Fold1(const rna_t& rna, fold_state_t* fold_state) {
  return FoldInternal(rna, fold_state, ComputeTables1);
}

folded_rna_t Fold2(const rna_t& rna, fold_state_t* fold_state) {
  return FoldInternal(rna, fold_state, ComputeTables2);
}

folded_rna_t Fold3(const rna_t& rna, fold_state_t* fold_state) {
  return FoldInternal(rna, fold_state, ComputeTables3);
}

energy_t FastTwoLoop(int ost, int oen, int ist, int ien) {
  int toplen = ist - ost - 1, botlen = oen - ien - 1;
  if (toplen == 0 && botlen == 0)
    return g_stack[r[ost]][r[ist]][r[ien]][r[oen]];
  if (toplen == 0 || botlen == 0)
    return energy::Bulge(ost, oen, ist, ien);
  if (toplen == 1 && botlen == 1)
    return g_internal_1x1[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[oen]];
  if (toplen == 1 && botlen == 2)
    return g_internal_1x2[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]][r[oen]];
  if (toplen == 2 && botlen == 1)
    return g_internal_1x2[r[ien]][r[ien + 1]][r[oen]][r[ost]][r[ost + 1]][r[ost + 2]][r[ist]];
  if (toplen == 2 && botlen == 2)
    return g_internal_2x2[r[ost]][r[ost + 1]][r[ost + 2]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]][r[oen]];

  static_assert(constants::TWOLOOP_MAX_SZ <= INITIATION_CACHE_SZ, "initiation cache not large enough");
  energy_t energy =
      g_internal_init[toplen + botlen] +
          std::min(std::abs(toplen - botlen) * g_internal_asym, constants::NINIO_MAX_ASYM);

  energy += energy::InternalLoopAuGuPenalty(r[ost], r[oen]);
  energy += energy::InternalLoopAuGuPenalty(r[ist], r[ien]);

  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2))
    energy += g_internal_2x3_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
        g_internal_2x3_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];
  else if (toplen != 1 && botlen != 1)
    energy += g_internal_other_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
        g_internal_other_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];

  return energy;
}

std::vector<hairpin_precomp_t> PrecomputeFastHairpin() {
  assert(r.size() > 0);
  std::vector<hairpin_precomp_t> precompute(r.size());
  std::string rna_str = parsing::RnaToString(r);
  for (const auto& hairpinpair : g_hairpin_e) {
    const auto& str = hairpinpair.first;
    verify_expr(str.size() - 2 <= MAX_SPECIAL_HAIRPIN_SZ, "need to increase MAX_SPECIAL_HAIRPIN_SZ");
    auto pos = rna_str.find(str, 0);
    while (pos != std::string::npos) {
      precompute[pos].special[str.size() - 2] = hairpinpair.second;
      pos = rna_str.find(str, pos + 1);
    }
  }
  int N = int(r.size());
  precompute[N - 1].num_c = int(r[N - 1] == C);
  for (int i = N - 2; i >= 0; --i) {
    if (r[i] == C)
      precompute[i].num_c = precompute[i + 1].num_c + 1;
  }
  return precompute;
}

energy_t FastHairpin(int st, int en, const std::vector<hairpin_precomp_t>& precomp) {
  int length = en - st - 1;
  assert(length >= constants::HAIRPIN_MIN_SZ);
  if (length <= MAX_SPECIAL_HAIRPIN_SZ && precomp[st].special[length] != MAX_E)
    return precomp[st].special[length];
  base_t stb = r[st], st1b = r[st + 1], en1b = r[en - 1], enb = r[en];
  energy_t energy = energy::HairpinInitiation(length) + energy::AuGuPenalty(stb, enb);

  bool all_c = precomp[st + 1].num_c >= length;

  if (length == 3) {
    if (all_c)
      energy += g_hairpin_c3_loop;
    return energy;
  }
  energy += g_terminal[r[st]][st1b][en1b][r[en]];

  if ((st1b == U && en1b == U) || (st1b == G && en1b == A))
    energy += g_hairpin_uu_ga_first_mismatch;
  if (st1b == G && en1b == G)
    energy += g_hairpin_gg_first_mismatch;
  if (all_c)
    energy += g_hairpin_all_c_a * length + g_hairpin_all_c_b;
  if (stb == G && enb == U && st >= 2 && r[st - 1] == G && r[st - 2] == G)
    energy += g_hairpin_special_gu_closure;

  return energy;
}

energy_t MinEnergy(const energy_t* energy, std::size_t size) {
  energy_t min = energy[0];
  for (int i = 0; i < int(size / sizeof(energy_t)); ++i)
    min = std::min(min, energy[i]);
  return min;
}

int MaxNumContiguous(const rna_t& rna) {
  energy_t num_contig = 0;
  energy_t max_num_contig = 0;
  base_t prev = -1;
  for (auto b : rna) {
    if (b == prev) num_contig++;
    else num_contig = 1;
    prev = b;
    max_num_contig = std::max(max_num_contig, num_contig);
  }
  return max_num_contig;
}

void InitFold() {
  // Initialise fast AUGU branch table
  for (base_t i = 0; i < 4; ++i)
    for (base_t j = 0; j < 4; ++j)
      g_augubranch[i][j] = g_multiloop_hack_b + energy::AuGuPenalty(i, j);

  auto min_stack = MinEnergy(&g_stack[0][0][0][0], sizeof(g_stack));

  // Non continuous (-2.1), -4 for WC, -16 for terminal mismatch.
  g_min_mismatch_coax = g_coax_mismatch_non_contiguous +
      std::min(g_coax_mismatch_gu_bonus, g_coax_mismatch_wc_bonus) +
      MinEnergy(&g_terminal[0][0][0][0], sizeof(g_terminal));
  // Minimum of all stacking params.
  g_min_flush_coax = min_stack;

  energy_t min_internal = MinEnergy(&g_internal_1x1[0][0][0][0][0][0], sizeof(g_internal_1x1));
  min_internal = std::min(min_internal,
      MinEnergy(&g_internal_1x2[0][0][0][0][0][0][0], sizeof(g_internal_1x2)));
  min_internal = std::min(min_internal,
      MinEnergy(&g_internal_2x2[0][0][0][0][0][0][0][0], sizeof(g_internal_2x2)));
  verify_expr(g_internal_asym >= 0,
      "min_internal optimisation does not work for negative asymmetry penalties");
  auto min_mismatch = 2 * std::min(
      MinEnergy(&g_internal_2x3_mismatch[0][0][0][0], sizeof(g_internal_2x3_mismatch)),
      MinEnergy(&g_internal_other_mismatch[0][0][0][0], sizeof(g_internal_other_mismatch)));
  auto min_internal_init = MinEnergy(&g_internal_init[4], sizeof(g_internal_init) - 4 * sizeof(g_internal_init[0]));
  min_internal = std::min(min_internal,
      min_internal_init + std::min(2 * g_internal_augu_penalty, 0) + min_mismatch);

  auto min_bulge_init = MinEnergy(&g_bulge_init[1], sizeof(g_bulge_init) - sizeof(g_bulge_init[0]));

  energy_t states_bonus = -energy_t(round(10.0 * constants::R * constants::T * log(MaxNumContiguous(r))));
  energy_t min_bulge = min_bulge_init + std::min(2 * g_augu_penalty, 0) +
      min_stack + std::min(g_bulge_special_c, 0) + states_bonus;
  g_min_twoloop_not_stack = std::min(min_bulge, min_internal);
}


fold::fold_fn_t* FoldFunctionFromArgParse(const ArgParse& argparse) {
  fold::fold_fn_t* fold_fn = nullptr;
  auto opt = argparse.GetOption("alg");
  if (opt == "0")
    fold_fn = &fold::Fold0;
  else if (opt == "1")
    fold_fn = &fold::Fold1;
  else if (opt == "2")
    fold_fn = &fold::Fold2;
  else if (opt == "3")
    fold_fn = &fold::Fold3;
  else if (opt == "brute")
    fold_fn = &fold::FoldBruteForce;
  else
    verify_expr(false, "unknown fold option");
  return fold_fn;
}

}
}
