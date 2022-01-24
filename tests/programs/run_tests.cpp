// Copyright 2016 Eliot Courtney.
#include "common_test.h"
#include "gtest/gtest.h"

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  mrna::ArgParse args(mrna::energy::ENERGY_OPTS);
  args.ParseOrExit(argc, argv);
  mrna::g_em = mrna::energy::EnergyModel::FromArgParse(args);
  return RUN_ALL_TESTS();
}
