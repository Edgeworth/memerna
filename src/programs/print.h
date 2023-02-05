// Copyright 2022 Eliot Courtney.
#ifndef PROGRAMS_PRINT_H_
#define PROGRAMS_PRINT_H_

#include "compute/partition/partition.h"
#include "models/t04/part/part.h"

void PrintBoltzProbs(const mrna::BoltzProbs& p);

void PrintPartition(const mrna::part::Part& p);

#endif  // PROGRAMS_PRINT_H_
