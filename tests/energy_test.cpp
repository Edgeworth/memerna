#include <energy.h>
#include <globals.h>
#include "gtest/gtest.h"

using namespace memerna;

TEST(HairpinEnergy, SixLength) {
  r = {C, A, C, A, A, A, A, A, A, A, U, G, U, G};
  EXPECT_EQ(46, energy::HairpinEnergy(3, 10));
}
