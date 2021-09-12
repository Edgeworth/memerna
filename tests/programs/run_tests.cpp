// Copyright 2016 E.
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
