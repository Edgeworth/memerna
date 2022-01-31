// Copyright 2022 E.
#include "compute/energy/boltzmann_model.h"

#include <string>
#include <type_traits>
#include <utility>

namespace mrna::energy {

namespace {

template <typename T>
constexpr auto Decay(T& a) {
  return reinterpret_cast<std::remove_all_extents_t<T>*>(&a);
}

void FillBoltzArray(BoltzEnergy* output, const Energy* input, int elements) {
  for (int i = 0; i < elements; ++i) output[i] = Boltz(input[i]);
}

}  // namespace

BoltzEnergyModel::BoltzEnergyModel(EnergyModelPtr em) : em_(em) {
#define FILL_BOLTZ(name)                                                                  \
  static_assert(/* NOLINTNEXTLINE */                                                      \
      sizeof(name) / sizeof(*Decay(name)) == sizeof(em->name) / sizeof(*Decay(em->name)), \
      "BoltzEnergyModel does not match EnergyModel");                                     \
  /* NOLINTNEXTLINE */                                                                    \
  FillBoltzArray(Decay(name), Decay(em->name), sizeof(name) / sizeof(*Decay(name)));

  FILL_BOLTZ(stack);
  FILL_BOLTZ(terminal);
  FILL_BOLTZ(internal_init);
  FILL_BOLTZ(internal_1x1);
  FILL_BOLTZ(internal_1x2);
  FILL_BOLTZ(internal_2x2);
  FILL_BOLTZ(internal_2x3_mismatch);
  FILL_BOLTZ(internal_other_mismatch);
  FILL_BOLTZ(internal_asym);
  FILL_BOLTZ(internal_augu_penalty);
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
  FILL_BOLTZ(augu_penalty);

  for (const auto& kv : em->hairpin) hairpin[kv.first] = Boltz(kv.second);
#undef FILL_BOLTZ
}

}  // namespace mrna::energy
