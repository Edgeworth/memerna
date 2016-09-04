#include "gtest/gtest.h"
#include "base.h"
#include "energy/energy_model.h"
#include "ctds_test.h"

using namespace memerna;

int main(int argc, char** argv) {
  energy::InitCtdsTest();
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
