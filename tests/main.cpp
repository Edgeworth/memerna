#include "gtest/gtest.h"
#include "base.h"

using namespace memerna;

int main(int argc, char** argv) {
  LoadEnergyModelFromDataDir();
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
