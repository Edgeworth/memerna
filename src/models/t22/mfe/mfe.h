// Copyright 2023 E.

#ifndef MODELS_T22_MFE_MFE_H_
#define MODELS_T22_MFE_MFE_H_

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

void MfeSlowest(const Primary& r, const Model::Ptr& em, DpState& state);

inline Energy MfeExterior(const Primary& r, const Model::Ptr& em, DpState& state) {
  return t04::MfeExterior(r, *em, state.t04);
}

}  // namespace mrna::md::t22

#endif  // MODELS_T22_MFE_MFE_H_
