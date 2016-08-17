#include <stack>
#include "parsing.h"
#include "fold/fold.h"
#include "fold/slow_fold.h"
#include "fold/traceback.h"
#include "fold/fold1.h"
#include "fold/fold2.h"
#include "fold/fold_globals.h"

namespace memerna {
namespace fold {

using constants::MAX_E;

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

  if (IsAuGu(r[ost], r[oen]))
    energy += g_internal_augu_penalty;
  if (IsAuGu(r[ist], r[ien]))
    energy += g_internal_augu_penalty;

  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2))
    energy += g_internal_2x3_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
        g_internal_2x3_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];
  else if (toplen != 1 && botlen != 1)
    energy += g_internal_other_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
        g_internal_other_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];

  return energy;
}

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

std::vector<energy_t[MAX_SPECIAL_HAIRPIN_SZ + 1]> PrecomputeFastHairpin() {
  std::vector<energy_t[MAX_SPECIAL_HAIRPIN_SZ + 1]> precompute(r.size());
  memset(precompute.data(), MAX_E & 0xFF, precompute.size() * sizeof(precompute[0]));
  std::string rna_str = parsing::RnaToString(r);
  for (const auto& hairpinpair : g_hairpin_e) {
    const auto& str = hairpinpair.first;
    auto pos = rna_str.find(str, 0);
    while (pos != std::string::npos) {
      precompute[pos][str.size()] = hairpinpair.second;
      pos = rna_str.find(str, pos+1);
    }
  }
  return precompute;
}

energy_t FastHairpin(int st, int en) {
  // Subtract two for the initiating base pair.
  int length = en - st - 1;
  if (length < 3) return constants::MAX_E;  // Disallowed by T04.
  energy_t energy = energy::HairpinInitiation(length);
  // Apply AU penalty if necessary (N.B. not for special hairpin sequences).
  if (IsAuGu(r[st], r[en])) {
    energy += g_augu_penalty;
  }

  // T04 says hairpin loops with all C bases inside them are treated specially.
  bool all_c = true;
  for (int i = st + 1; i <= en - 1; ++i) {
    if (r[i] != C) all_c = false;
  }

  if (length == 3) {
    if (all_c) {
      energy += g_hairpin_c3_loop;
    }
    return energy;
  }
  base_t left = r[st + 1], right = r[en - 1];
  energy += g_terminal[r[st]][left][right][r[en]];
  if ((left == U && right == U) || (left == G_b && right == A_b))
    energy += g_hairpin_uu_ga_first_mismatch;
  if (IsPairOf(left, right, G_b, G_b))
    energy += g_hairpin_gg_first_mismatch;
  if (all_c)
    energy += g_hairpin_all_c_a * length + g_hairpin_all_c_b;
  if (r[st] == G && r[en] == U && st >= 2 && r[st - 1] == G && r[st - 2] == G)
    energy += g_hairpin_special_gu_closure;

  return energy;
}

void InitFold() {
  if (g_fold_init)
    return;
  // Initialise fast AUGU branch table
  for (base_t i = 0; i < 4; ++i)
    for (base_t j = 0; j < 4; ++j)
      g_augubranch[i][j] = g_multiloop_hack_b + energy::AuGuPenalty(i, j);

  g_fold_init = true;
}

}
}
