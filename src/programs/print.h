// Copyright 2022 Eliot Courtney.
#ifndef PROGRAMS_PRINT_H_
#define PROGRAMS_PRINT_H_

#include "model/part.h"

void PrintBoltzProbs(const mrna::BoltzProbs& p);

void PrintPartition(const mrna::BoltzSums& p);

#endif  // PROGRAMS_PRINT_H_
