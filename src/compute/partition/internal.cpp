// Copyright 2016 E.
#include <algorithm>

#include "compute/energy/fast_energy.h"
#include "compute/energy/globals.h"
#include "compute/partition/globals.h"
#include "compute/partition/partition.h"

namespace mrna::partition::internal {

namespace {

template <typename T>
constexpr auto Decay(T& a) {
  return reinterpret_cast<std::remove_all_extents_t<T>*>(&a);
}

void FillPenergyArray(penergy_t* output, const energy_t* input, int elements) {
  for (int i = 0; i < elements; ++i) output[i] = energy::Boltzmann(input[i]);
}

}  // namespace

penergy_model_t::penergy_model_t(const energy::EnergyModel& em) {
#define FILL_PENERGY(name)                                                                        \
  static_assert(sizeof(name) / sizeof(*Decay(name)) == sizeof(em.name) / sizeof(*Decay(em.name)), \
      "penergy_model_t does not match EnergyModel");                                              \
  FillPenergyArray(Decay(name), Decay(em.name), sizeof(name) / sizeof(*Decay(name)));

  FILL_PENERGY(stack);
  FILL_PENERGY(terminal);
  FILL_PENERGY(internal_init);
  FILL_PENERGY(internal_1x1);
  FILL_PENERGY(internal_1x2);
  FILL_PENERGY(internal_2x2);
  FILL_PENERGY(internal_2x3_mismatch);
  FILL_PENERGY(internal_other_mismatch);
  FILL_PENERGY(internal_asym);
  FILL_PENERGY(internal_augu_penalty);
  FILL_PENERGY(bulge_init);
  FILL_PENERGY(bulge_special_c);
  FILL_PENERGY(hairpin_init);
  FILL_PENERGY(hairpin_uu_ga_first_mismatch);
  FILL_PENERGY(hairpin_gg_first_mismatch);
  FILL_PENERGY(hairpin_special_gu_closure);
  FILL_PENERGY(hairpin_c3_loop);
  FILL_PENERGY(hairpin_all_c_a);
  FILL_PENERGY(hairpin_all_c_b);
  FILL_PENERGY(multiloop_hack_a);
  FILL_PENERGY(multiloop_hack_b);
  FILL_PENERGY(dangle5);
  FILL_PENERGY(dangle3);
  FILL_PENERGY(coax_mismatch_non_contiguous);
  FILL_PENERGY(coax_mismatch_wc_bonus);
  FILL_PENERGY(coax_mismatch_gu_bonus);
  FILL_PENERGY(augu_penalty);

  for (const auto& kv : em.hairpin) hairpin[kv.first] = energy::Boltzmann(kv.second);
#undef FILL_PENERGY
}

penergy_t FastHairpin(int st, int en) {
  int length = en - st - 1;
  assert(length >= HAIRPIN_MIN_SZ);
  if (length <= energy::MAX_SPECIAL_HAIRPIN_SZ && gppc.hairpin[st].special[length] > -1)
    return gppc.hairpin[st].special[length];
  base_t stb = gr[st], st1b = gr[st + 1], en1b = gr[en - 1], enb = gr[en];
  penergy_t energy =
      energy::Boltzmann(energy::gem.HairpinInitiation(length) + energy::gem.AuGuPenalty(stb, enb));

  bool all_c = gppc.hairpin[st + 1].num_c >= length;

  if (length == 3) {
    if (all_c) energy *= gppc.em.hairpin_c3_loop;
    return energy;
  }
  energy *= gppc.em.terminal[gr[st]][st1b][en1b][gr[en]];

  if ((st1b == U && en1b == U) || (st1b == G && en1b == A))
    energy *= gppc.em.hairpin_uu_ga_first_mismatch;
  if (st1b == G && en1b == G) energy *= gppc.em.hairpin_gg_first_mismatch;
  if (all_c)
    energy *= energy::Boltzmann(energy::gem.hairpin_all_c_a * length) * gppc.em.hairpin_all_c_b;
  if (stb == G && enb == U && st >= 2 && gr[st - 1] == G && gr[st - 2] == G)
    energy *= gppc.em.hairpin_special_gu_closure;

  return energy;
}

penergy_t FastTwoLoop(int ost, int oen, int ist, int ien) {
  int toplen = ist - ost - 1, botlen = oen - ien - 1;
  if (toplen == 0 && botlen == 0) return gppc.em.stack[gr[ost]][gr[ist]][gr[ien]][gr[oen]];
  if (toplen == 0 || botlen == 0)
    return energy::Boltzmann(energy::gem.Bulge(gr, ost, oen, ist, ien));
  if (toplen == 1 && botlen == 1)
    return gppc.em.internal_1x1[gr[ost]][gr[ost + 1]][gr[ist]][gr[ien]][gr[ien + 1]][gr[oen]];
  if (toplen == 1 && botlen == 2)
    return gppc.em
        .internal_1x2[gr[ost]][gr[ost + 1]][gr[ist]][gr[ien]][gr[ien + 1]][gr[ien + 2]][gr[oen]];
  if (toplen == 2 && botlen == 1)
    return gppc.em
        .internal_1x2[gr[ien]][gr[ien + 1]][gr[oen]][gr[ost]][gr[ost + 1]][gr[ost + 2]][gr[ist]];
  if (toplen == 2 && botlen == 2)
    return gppc.em.internal_2x2[gr[ost]][gr[ost + 1]][gr[ost + 2]][gr[ist]][gr[ien]][gr[ien + 1]]
                               [gr[ien + 2]][gr[oen]];

  static_assert(TWOLOOP_MAX_SZ <= energy::EnergyModel::INITIATION_CACHE_SZ,
      "initiation cache not large enough");
  assert(toplen + botlen < energy::EnergyModel::INITIATION_CACHE_SZ);
  penergy_t energy = gppc.em.internal_init[toplen + botlen] *
      energy::Boltzmann(
          std::min(std::abs(toplen - botlen) * energy::gem.internal_asym, NINIO_MAX_ASYM));

  energy *= gppc.em.InternalLoopAuGuPenalty(gr[ost], gr[oen]);
  energy *= gppc.em.InternalLoopAuGuPenalty(gr[ist], gr[ien]);

  if ((toplen == 2 && botlen == 3) || (toplen == 3 && botlen == 2))
    energy *= gppc.em.internal_2x3_mismatch[gr[ost]][gr[ost + 1]][gr[oen - 1]][gr[oen]] *
        gppc.em.internal_2x3_mismatch[gr[ien]][gr[ien + 1]][gr[ist - 1]][gr[ist]];
  else if (toplen != 1 && botlen != 1)
    energy *= gppc.em.internal_other_mismatch[gr[ost]][gr[ost + 1]][gr[oen - 1]][gr[oen]] *
        gppc.em.internal_other_mismatch[gr[ien]][gr[ien + 1]][gr[ist - 1]][gr[ist]];

  return energy;
}

precomp_t PrecomputeData(const primary_t& r, const energy::EnergyModel& em) {
  precomp_t pc;
  pc.em = penergy_model_t(em);
  for (base_t i = 0; i < 4; ++i)
    for (base_t j = 0; j < 4; ++j)
      pc.augubranch[i][j] = pc.em.multiloop_hack_b * pc.em.AuGuPenalty(i, j);
  pc.hairpin = energy::PrecomputeHairpin<energy::hairpin_precomp_t<penergy_t, -1>>(r, pc.em);
  return pc;
}

}  // namespace mrna::partition::internal
