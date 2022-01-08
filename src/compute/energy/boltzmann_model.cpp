// Copyright 2022 E.
#include "compute/energy/boltzmann_model.h"

#include <algorithm>
#include <cassert>
#include <memory>
#include <string>

#include "compute/energy/structure.h"
#include "model/parsing.h"
#include "util/string.h"

namespace mrna::energy {

namespace {

template <typename T>
constexpr auto Decay(T& a) {
  return reinterpret_cast<std::remove_all_extents_t<T>*>(&a);
}

void FillPenergyArray(BoltzEnergy* output, const Energy* input, int elements) {
  for (int i = 0; i < elements; ++i) output[i] = Boltz(input[i]);
}

}  // namespace

BoltzEnergyModel::BoltzEnergyModel(EnergyModel em) {
#define FILL_PENERGY(name)                                                                        \
  static_assert(sizeof(name) / sizeof(*Decay(name)) == sizeof(em.name) / sizeof(*Decay(em.name)), \
      "BoltzEnergyModel does not match EnergyModel");                                             \
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

  for (const auto& kv : em.hairpin) hairpin[kv.first] = Boltz(kv.second);
#undef FILL_PENERGY
}

}  // namespace mrna::energy
