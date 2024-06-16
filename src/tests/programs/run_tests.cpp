// Copyright 2016 Eliot Courtney.
#include <gtest/gtest.h>

#include "tests/init.h"
#include "util/argparse.h"
#include "util/error.h"

int main(int argc, char** argv) {
  mrna::InitProgram();
  testing::InitGoogleTest(&argc, argv);
  mrna::ArgParse args;
  mrna::RegisterOptsBackendCfg(&args);
  args.ParseOrExit(argc, argv);
  mrna::InitTest(args.Get(mrna::OPT_MEMERNA_DATA));

  return RUN_ALL_TESTS();
}
