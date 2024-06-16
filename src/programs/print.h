// Copyright 2022 Eliot Courtney.
#ifndef PROGRAMS_PRINT_H_
#define PROGRAMS_PRINT_H_

#include "model/pfn.h"

void PrintBoltzProbs(const mrna::BoltzProbs& p);

void PrintPfn(const mrna::BoltzSums& p);

#endif  // PROGRAMS_PRINT_H_
