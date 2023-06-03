// Copyright 2016 Eliot Courtney.
#include <gtest/gtest.h>

#include <ios>

#include "api/energy/energy_cfg.h"
#include "common_test.h"
#include "util/argparse.h"

int main(int argc, char** argv) {
  std::ios_base::sync_with_stdio(false);
  testing::InitGoogleTest(&argc, argv);
  mrna::ArgParse args;
  mrna::erg::RegisterOptsEnergyModel(&args);
  args.ParseOrExit(argc, argv);
  mrna::InitTest(args.Get(mrna::erg::OPT_MEMERNA_DATA));

  return RUN_ALL_TESTS();
}
