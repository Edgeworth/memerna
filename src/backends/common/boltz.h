// Copyright 2023 Eliot Courtney.
#ifndef BACKENDS_COMMON_BOLTZ_H_
#define BACKENDS_COMMON_BOLTZ_H_

#include <cstddef>
#include <type_traits>

#include "model/energy.h"
#include "util/util.h"  // IWYU pragma: keep

namespace mrna::md {

template <typename T, typename U>
struct SameArrayShape : std::bool_constant<!std::is_array_v<T> && !std::is_array_v<U>> {};

template <typename T, typename U, std::size_t N, std::size_t M>
struct SameArrayShape<T[N], U[M]> : std::bool_constant<N == M && SameArrayShape<T, U>::value> {};

inline void FillBoltzArray(BoltzEnergy* output, const Energy* input, std::size_t elements) {
  for (std::size_t i = 0; i < elements; ++i) output[i] = Boltz(input[i]);
}

template <typename Output, typename Input>
inline void FillBoltzArray(Output& output, const Input& input) {
  static_assert(SameArrayShape<Output, Input>::value, "BoltzModel does not match Model");
  auto* output_data = Decay(output);
  const auto* input_data = Decay(input);
  FillBoltzArray(output_data, input_data, sizeof(output) / sizeof(*output_data));
}

#define FILL_BOLTZ(bm, name) FillBoltzArray((bm).name, (bm).m().name)

}  // namespace mrna::md

#endif  // BACKENDS_COMMON_BOLTZ_H_
