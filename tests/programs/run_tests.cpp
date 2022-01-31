// Copyright 2016 Eliot Courtney.
#include <cassert>

#include "common_test.h"
#include "compute/energy/config.h"
#include "compute/energy/model.h"
#include "gtest/gtest.h"
#include "util/argparse.h"

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  mrna::ArgParse args;
  mrna::energy::RegisterOpts(&args);
  args.ParseOrExit(argc, argv);

  mrna::t04 = mrna::energy::EnergyModel::FromArgParse(args);
  assert(mrna::t04->Checksum() == mrna::T04_MODEL_HASH);

  mrna::g_em[0] = mrna::t04;
  for (int i = 1; i < mrna::NUM_TEST_MODELS; ++i)
    mrna::g_em[i] = mrna::energy::EnergyModel::Random(i);

  return RUN_ALL_TESTS();
}
