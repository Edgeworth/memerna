// Copyright 2023 Eliot Courtney.

#ifndef COMPUTE_MFE_T22_MFE_H_
#define COMPUTE_MFE_T22_MFE_H_

#include "compute/energy/t22/model.h"
#include "compute/mfe/t04/dp.h"
#include "compute/mfe/t04/mfe.h"
#include "compute/mfe/t22/dp.h"
#include "model/constants.h"
#include "model/primary.h"

namespace mrna::md::t22::mfe {

struct DpState {
  // T04 state is reused.
  t04::DpState t04;
  // DP for stack length.
  Array3D<Energy> stack;
};

void MfeSlowest(const Primary& r, const erg::t22::Model::Ptr& em, DpState& state);

inline void MfeExterior(const Primary& r, const erg::t22::Model::Ptr& em, DpState& state) {
  t04::MfeExterior(r, *em, state.t04);
}

}  // namespace mrna::md::t22::mfe

#endif  // COMPUTE_MFE_T22_MFE_H_
