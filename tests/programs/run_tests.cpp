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
#include "common_test.h"
#include "energy/load_model.h"
#include "gtest/gtest.h"

using namespace memerna;

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  ArgParse argparse(energy::ENERGY_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  g_em = energy::LoadEnergyModelFromArgParse(argparse);
  g_ems.push_back(g_em);
  for (int_fast32_t i = 0; i < 4; ++i) g_ems.push_back(energy::LoadRandomEnergyModel(i));
  return RUN_ALL_TESTS();
}
