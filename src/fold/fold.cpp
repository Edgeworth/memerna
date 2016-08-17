#include <stack>
#include "fold/fold.h"
#include "fold/slow_fold.h"
#include "fold/traceback.h"
#include "fold/fold1.h"
#include "fold/fold2.h"

namespace memerna {
namespace fold {

energy_t FastTwoLoop(int ost, int oen, int ist, int ien) {
  int toplen = ist - ost - 1, botlen = oen - ien - 1;
  if (toplen == 0 && botlen == 0)
    return stacking_e[r[ost]][r[ist]][r[ien]][r[oen]];
  if (toplen == 0 || botlen == 0)
    return energy::BulgeEnergy(ost, oen, ist, ien);
  if (toplen == 1 && botlen == 1)
    return internal_1x1[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[oen]];
  if (toplen == 1 && botlen == 2)
    return internal_1x2[r[ost]][r[ost + 1]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]][r[oen]];
  if (toplen == 2 && botlen == 1)
    return internal_1x2[r[ien]][r[ien + 1]][r[oen]][r[ost]][r[ost + 1]][r[ost + 2]][r[ist]];
  if (toplen == 2 && botlen == 2)
    return internal_2x2[r[ost]][r[ost + 1]][r[ost + 2]][r[ist]][r[ien]][r[ien + 1]][r[ien + 2]][r[oen]];

  static_assert(constants::TWOLOOP_MAX_SZ <= INITIATION_CACHE_SZ, "initiation cache not large enough");
  energy_t energy = internal_init[toplen + botlen] + std::min(std::abs(toplen - botlen) * internal_asym, constants::NINIO_MAX_ASYM);

  if (IsAuGu(r[ost], r[oen]))
    energy += internal_augu_penalty;
  if (IsAuGu(r[ist], r[ien]))
    energy += internal_augu_penalty;

  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2))
    energy += internal_2x3_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
        internal_2x3_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];
  else if (toplen != 1 && botlen != 1)
    energy += internal_other_mismatch[r[ost]][r[ost + 1]][r[oen - 1]][r[oen]] +
        internal_other_mismatch[r[ien]][r[ien + 1]][r[ist - 1]][r[ist]];

  return energy;
}

template<typename T>
energy_t FoldInternal(T ComputeTables) {
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

}
}
