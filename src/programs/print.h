// Copyright 2022 Eliot Courtney.
#ifndef PROGRAMS_PRINT_H_
#define PROGRAMS_PRINT_H_

#include <span>
#include <string>

#include "api/ctx/backend.h"
#include "api/energy/pseudofree_cfg.h"
#include "model/energy.h"
#include "model/pfn.h"
#include "model/primary.h"

namespace mrna {

void PrintBoltzProbs(const BoltzProbs& p);

void PrintPfn(const BoltzSums& p);

void PrintExtraPfnTables(const Primary& r, const PfnTables& pfn, const BackendModelPtr& model,
    const erg::PseudofreeCfg& pf);

std::string FormatPseudofreeEnergies(std::span<const Energy> energies);

}  // namespace mrna

#endif  // PROGRAMS_PRINT_H_
