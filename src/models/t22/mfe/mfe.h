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
  // DP for holding the best answer given (st, en) is paired but we can't
  // start a stack of length 2 or more from here.
  Array2D<Energy> nostack;
  // DP for computing the best stack of a given length.
  Array3D<Energy> penult;
};

// Use int16_t here to save memory.
struct PenultimateIndex {
  Index st, en, len;

  PenultimateIndex(int st_, int en_, int len_) : st(Index(st_)), en(Index(en_)), len(Index(len_)) {
    assert(st_ == st && en_ == en && len == len_);
  }

  constexpr auto operator<=>(const PenultimateIndex&) const = default;
};

struct NoStackIndex {
  Index st, en;

  NoStackIndex(int st_, int en_) : st(Index(st_)), en(Index(en_)) {
    assert(st_ == st && en_ == en);
  }

  constexpr auto operator<=>(const NoStackIndex&) const = default;
};

// Index into the DP tables.
using DpIndex = std::variant<t04::DpIndex, NoStackIndex, PenultimateIndex>;

void MfeDebug(const Primary& r, const Model::Ptr& em, DpState& state);

Energy MfeExterior(const Primary& r, const Model::Ptr& em, DpState& state);

}  // namespace mrna::md::t22

#endif  // MODELS_T22_MFE_MFE_H_
