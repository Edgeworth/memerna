#include "gtest/gtest.h"
#include "base.h"
#include "energy/energy_model.h"
#include "common_test.h"
#include "ctds_test.h"

using namespace memerna;

int main(int argc, char** argv) {
  energy::LoadEnergyModelFromDataDir(ENERGY_MODEL_PATH);
  energy::InitCtdsTest();
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
