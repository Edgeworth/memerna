// Copyright 2022 Eliot Courtney.
#ifndef COMPUTE_BOLTZ_DP_H_
#define COMPUTE_BOLTZ_DP_H_

#include "model/base.h"
#include "model/ctd.h"
#include "util/array.h"

namespace mrna {

// DP arrays
enum : int8_t { PT_P, PT_U, PT_U2, PT_U_WC, PT_U_GU, PT_U_RC, PT_SIZE };

using BoltzDpArray = Array2D1S<BoltzEnergy, PT_SIZE>;

enum : int8_t {
  PTEXT_R,
  PTEXT_L,
  PTEXT_R_WC,  // Must start with a branch not involved in an interaction that is Watson-Crick
  PTEXT_R_GU,  // Must start with a branch not involved in an interaction that is GU
  // Must start with a branch, that branch is involved backwards in a right coaxial stack.
  PTEXT_R_RC,
  PTEXT_L_WC,
  PTEXT_L_GU,
  PTEXT_L_LCOAX,
  PTEXT_SIZE
};

using BoltzExtArray = Array1D1S<BoltzEnergy, PTEXT_SIZE>;

// Holds the Boltzmann sums used in the partition function.ss
using BoltzSums = Array2D<BoltzEnergy>;

using BoltzProbs = Array2D<BoltzEnergy>;

}  // namespace mrna

#endif  // COMPUTE_BOLTZ_DP_H_
