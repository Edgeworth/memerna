// Copyright 2016, E.
//
// This file is part of memerna.
//
// memerna is free software: you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// memerna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with memerna.
// If not, see <http://www.gnu.org/licenses/>.
#ifndef MEMERNA_MEMERNA_H
#define MEMERNA_MEMERNA_H

#include "bridge/bridge.h"
#include "common.h"
#include "context.h"

namespace memerna {
namespace bridge {

// Note that only one energy model can be loaded at a time.
class Memerna : public RnaPackage {
public:
  Memerna(energy::EnergyModelPtr em_, context_opt_t options_)
      : em(em_), options(options_) {}

  Memerna(const Memerna&) = delete;
  Memerna& operator=(const Memerna&) = delete;

  virtual energy_t Efn(const secondary_t& secondary, std::string* desc = nullptr) const override;
  virtual computed_t Fold(const primary_t& r) const override;
  virtual int Suboptimal(fold::SuboptimalCallback fn,
      const primary_t& r, energy_t energy_delta) const override;
  virtual std::vector<computed_t> SuboptimalIntoVector(
      const primary_t& r, energy_t energy_delta) const override;
  virtual std::pair<partition::partition_t, partition::probabilities_t> Partition(
      const primary_t& r) const override;

private:
  const energy::EnergyModelPtr em;
  const context_opt_t options;
};
}
}

#endif  // MEMERNA_MEMERNA_H
