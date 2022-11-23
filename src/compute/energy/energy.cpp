// Copyright 2022 Eliot Courtney.
#include "compute/energy/energy.h"

#include <algorithm>
#include <utility>

#include "compute/energy/model.h"
#include "compute/energy/structure.h"

namespace mrna::erg {

int MaxNumContiguous(const Primary& r) {
  int num_contig = 0;
  int max_num_contig = 0;
  Base prev = -1;
  for (auto b : r) {
    if (b == prev)
      num_contig++;
    else
      num_contig = 1;
    prev = b;
    max_num_contig = std::max(max_num_contig, num_contig);
  }
  return max_num_contig;
}

EnergyResult::EnergyResult(Energy energy, Ctds ctd, std::unique_ptr<Structure> struc)
    : energy(energy), ctd(std::move(ctd)), struc(std::move(struc)) {}

EnergyResult::~EnergyResult() = default;

}  // namespace mrna::erg
