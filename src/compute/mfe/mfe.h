// Copyright 2016 E.
#ifndef COMPUTE_MFE_MFE_H_
#define COMPUTE_MFE_MFE_H_

#include "common.h"
#include "compute/mfe/globals.h"
#include "model/base.h"
#include "model/globals.h"

namespace mrna::mfe::internal {

struct cand_t {
  energy_t energy;
  int idx;
};

void ComputeTables0();
void ComputeTables1();
void ComputeTables2();
void ComputeTables3();
void ComputeExterior();
void Traceback();

}  // namespace mrna::mfe::internal

#endif  // COMPUTE_MFE_MFE_H_
