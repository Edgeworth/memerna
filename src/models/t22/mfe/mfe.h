// Copyright 2023 Eliot Courtney.

#ifndef MODELS_T22_MFE_MFE_H_
#define MODELS_T22_MFE_MFE_H_

#include "model/constants.h"
#include "model/primary.h"
#include "models/t04/mfe/dp.h"
#include "models/t04/mfe/mfe.h"
#include "models/t22/energy/model.h"
#include "models/t22/mfe/mfe.h"

namespace mrna::md::t22 {

struct DpState {
  // T04 state is reused.
  t04::DpState t04;
  // DP for stack length.
  Array3D<Energy> stack;
};

void MfeSlowest(const Primary& r, const Model::Ptr& em, DpState& state);

inline void MfeExterior(const Primary& r, const Model::Ptr& em, DpState& state) {
  t04::MfeExterior(r, *em, state.t04);
}

}  // namespace mrna::md::t22

#endif  // MODELS_T22_MFE_MFE_H_
