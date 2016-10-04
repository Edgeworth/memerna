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
#ifndef MEMERNA_RNASTRUCTURE_H
#define MEMERNA_RNASTRUCTURE_H

#include "bridge/bridge.h"
#include "common.h"

#include "miles_rnastructure/include/algorithm.h"
#include "miles_rnastructure/include/rna_library.h"

namespace memerna {
namespace bridge {

class Rnastructure : public RnaPackage {
public:
  Rnastructure(const std::string& data_path, bool use_lyngso_);
  Rnastructure(const Rnastructure&) = delete;
  Rnastructure& operator=(Rnastructure&) = delete;

  virtual energy_t Efn(const secondary_t& secondary, std::string* desc = nullptr) const override;
  virtual computed_t Fold(const primary_t& r) const override;
  virtual int Suboptimal(std::function<void(const computed_t&)> fn,
      const primary_t& r, energy_t energy_delta) const override;
  virtual std::vector<computed_t>  SuboptimalIntoVector(
      const primary_t& r, energy_t energy_delta) const override;

  computed_t FoldAndDpTable(const primary_t& r, dp_state_t* dp_state) const;

private:
  const std::unique_ptr<datatable> data;
  const bool use_lyngso;
};
}
}

#endif  // MEMERNA_RNASTRUCTURE_H
