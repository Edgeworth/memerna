// Copyright 2022 E.
#ifndef COMPUTE_PARTITION_T04_PARTITION_H_
#define COMPUTE_PARTITION_T04_PARTITION_H_

#include <cassert>
#include <tuple>

#include "model/constants.h"
#include "model/energy.h"
#include "model/primary.h"
#include "models/t04/energy/boltz_model.h"
#include "models/t04/part/part.h"
#include "util/array.h"

namespace mrna::md::t04::part {

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

struct PartState {
  BoltzDpArray dp;
  BoltzExtArray ext;
};

void PartitionSlowest(const Primary& r, const erg::t04::Model::Ptr& em, PartState& state);

void PartitionFastest(const Primary& r, const erg::t04::BoltzModel::Ptr& bem, PartState& state);

void PartitionExterior(const Primary& r, const erg::t04::Model& em, PartState& state);

}  // namespace mrna::md::t04::part

#endif  // COMPUTE_PARTITION_T04_PARTITION_H_