// Copyright 2016 E.
#include "common_test.h"
#include "energy/load_model.h"
#include "gtest/gtest.h"

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  mrna::ArgParse argparse(mrna::energy::ENERGY_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  mrna::g_em = mrna::energy::LoadEnergyModelFromArgParse(argparse);
  mrna::g_ems.push_back(mrna::g_em);
  for (int_fast32_t i = 0; i < 4; ++i)
    mrna::g_ems.push_back(mrna::energy::LoadRandomEnergyModel(i));
  return RUN_ALL_TESTS();
}
