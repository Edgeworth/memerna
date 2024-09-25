// Copyright 2022 Eliot Courtney.
#ifndef PROGRAMS_PRINT_H_
#define PROGRAMS_PRINT_H_

#include "backends/common/base/model_base.h"
#include "model/pfn.h"
#include "model/primary.h"

namespace mrna {

void PrintBoltzProbs(const BoltzProbs& p);

void PrintPfn(const BoltzSums& p);

void PrintInnerStackProbs(const Primary& r, const PfnTables& pfn, const md::base::ModelBase& model);

void PrintHelixProbs(const Primary& r, const PfnTables& pfn, const md::base::ModelBase& model);

}  // namespace mrna

#endif  // PROGRAMS_PRINT_H_
