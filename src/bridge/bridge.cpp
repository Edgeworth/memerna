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
#include "bridge/bridge.h"
#include "bridge/memerna.h"
#include "bridge/rnastructure.h"
#include "energy/load_model.h"

namespace memerna {
namespace bridge {

std::unique_ptr<RnaPackage> RnaPackageFromArgParse(const ArgParse& argparse) {
  verify_expr(argparse.HasFlag("r") + argparse.HasFlag("k") == 1,
      "require exactly one package flag\n%s", argparse.Usage().c_str());
  if (argparse.HasFlag("r")) {
    return std::unique_ptr<RnaPackage>(
        new Rnastructure("extern/miles_rnastructure/data_tables/", false));
  } else {
    return std::unique_ptr<RnaPackage>(new Memerna(
        energy::LoadEnergyModelFromArgParse(argparse), fold::ContextOptionsFromArgParse(argparse)));
  }
}
}
}
