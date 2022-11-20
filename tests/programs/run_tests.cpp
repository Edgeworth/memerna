// Copyright 2016 Eliot Courtney.
#include <cassert>
#include <memory>

#include "common_test.h"
#include "compute/energy/energy_cfg.h"
#include "compute/energy/model.h"
#include "gtest/gtest.h"
#include "util/argparse.h"
#include "util/error.h"

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  mrna::ArgParse args;
  mrna::energy::RegisterOpts(&args);
  args.ParseOrExit(argc, argv);

  // TODO(0): deal with precision stuff.
  mrna::t04 = mrna::energy::EnergyModel::FromArgParse(args);
  verify(mrna::t04->Checksum() == mrna::T04_P1_MODEL_HASH, "Expected T04 model hash %d, got %d",
      mrna::T04_P1_MODEL_HASH, mrna::t04->Checksum());

  mrna::g_em[0] = mrna::t04;
  for (int i = 1; i < mrna::NUM_TEST_MODELS; ++i)
    mrna::g_em[i] = mrna::energy::EnergyModel::Random(i);

  return RUN_ALL_TESTS();
}
