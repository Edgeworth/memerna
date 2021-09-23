// Copyright 2016 E.
#include "common_test.h"
#include "energy/load_model.h"
#include "gtest/gtest.h"

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  mrna::ArgParse argparse(mrna::energy::ENERGY_OPTIONS);
  argparse.ParseOrExit(argc, argv);
  mrna::g_em = mrna::energy::LoadEnergyModelFromArgParse(argparse);
  return RUN_ALL_TESTS();
}
