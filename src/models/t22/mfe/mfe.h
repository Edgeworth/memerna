// Copyright 2023 Eliot Courtney.

#ifndef MODELS_T22_MFE_MFE_H_
#define MODELS_T22_MFE_MFE_H_

#include <variant>

#include "model/constants.h"
#include "model/primary.h"
#include "models/t04/mfe/dp.h"
#include "models/t04/mfe/mfe.h"
#include "models/t22/energy/model.h"

namespace mrna::md::t22 {

struct DpState {
  // T04 state is reused.
  t04::DpState t04;
  // DP for penultimate stack length.
  Array3D<Energy> penult;
};

// Use int16_t here to save memory.
struct PenultimateIndex {
  int16_t st, en, len;

  PenultimateIndex(int st_, int en_, int len_)
      : st(int16_t(st_)), en(int16_t(en_)), len(int16_t(len_)) {
    assert(st_ == st && en_ == en && len == len_);
  }

  constexpr auto operator<=>(const PenultimateIndex&) const = default;
};

// Index into the DP tables.
using Index = std::variant<t04::Index, PenultimateIndex>;

void MfeSlowest(const Primary& r, const Model::Ptr& em, DpState& state);

inline Energy MfeExterior(const Primary& r, const Model::Ptr& em, DpState& state) {
  return t04::MfeExterior(r, *em, state.t04);
}

}  // namespace mrna::md::t22

#endif  // MODELS_T22_MFE_MFE_H_
