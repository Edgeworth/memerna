#include <stack>
#include "parsing.h"
#include "fold/fold.h"
#include "fold/slow_fold.h"
#include "fold/traceback.h"
#include "fold/fold1.h"
#include "fold/fold2.h"
#include "fold/fold3.h"
#include "fold/fold4.h"
#include "fold/fold_globals.h"

namespace memerna {
namespace fold {

using constants::MAX_E;

template<typename T>
energy_t FoldInternal(T ComputeTables) {
  InitFold();
  auto arr = ComputeTables();
  std::stack<std::tuple<int, int, int>> q;
  auto energy = TraceExterior(arr, q);
  TraceStructure(arr, q);
  return energy;
}

energy_t Fold() {
  return FoldInternal(ComputeTables2);
}

energy_t FoldSlow() {
  return FoldInternal(ComputeTablesSlow);
}

energy_t Fold1() {
  return FoldInternal(ComputeTables1);
}

energy_t Fold2() {
  return FoldInternal(ComputeTables2);
}

energy_t Fold3() {
  return FoldInternal(ComputeTables3);
}

energy_t Fold4() {
  return FoldInternal(ComputeTables4);
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

void InitFold() {
  // Initialise fast AUGU branch table
  for (base_t i = 0; i < 4; ++i)
    for (base_t j = 0; j < 4; ++j)
      g_augubranch[i][j] = g_multiloop_hack_b + energy::AuGuPenalty(i, j);
}

}
}
