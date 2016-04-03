#include <energy.h>
#include <globals.h>
#include <parsing.h>
#include "gtest/gtest.h"

using namespace memerna;

class EnergyTest : public testing::Test {
public:
    folded_rna_t kHairpin1 = parsing::ParseViennaRna("CACAAAAAAAUGUG", "((((......))))");
};


TEST_F(EnergyTest, SixLength) {
  r = kHairpin1.r;
  EXPECT_EQ(46, energy::HairpinEnergy(3, 10));
}

TEST_F(EnergyTest, NNDBExamples) {
  EXPECT_EQ(-14, energy::ComputeEnergy(kHairpin1));
}
