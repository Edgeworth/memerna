#include <energy.h>
#include <globals.h>
#include <parsing.h>
#include "gtest/gtest.h"

namespace memerna {
namespace energy {

class EnergyTest : public testing::Test {
public:
    folded_rna_t kHairpin1 = parsing::ParseViennaRna("CACAAAAAAAUGUG", "((((......))))");
};


TEST_F(EnergyTest, HairpinEnergy) {
  r = kHairpin1.r;
  EXPECT_EQ(AUGU_PENALTY + terminal_e[A][A][A][U] + HairpinInitiation(6), HairpinEnergy(3, 10));
}

TEST_F(EnergyTest, NNDBExamples) {
  EXPECT_EQ(
      stacking_e[C][A][U][G] + stacking_e[A][C][G][U] + stacking_e[C][A][U][G] + AUGU_PENALTY +
          terminal_e[A][A][A][U] + HairpinInitiation(6),
      ComputeEnergy(kHairpin1));
}

TEST_F(EnergyTest, BaseCases) {
  EXPECT_EQ(
      AUGU_PENALTY +  stacking_e[G][A][U][C] + hairpin_init[3],
      ComputeEnergy(parsing::ParseViennaRna("GAAAAUC", "((...))")));
  EXPECT_EQ(
      AUGU_PENALTY * 2 + stacking_e[G][A][U][U] + hairpin_init[3],
      ComputeEnergy(parsing::ParseViennaRna("GAAAAUU", "((...))")));
  EXPECT_EQ(AUGU_PENALTY * 2 + hairpin_init[3], ComputeEnergy(parsing::ParseViennaRna("AAAAAUA", ".(...).")));
  EXPECT_EQ(AUGU_PENALTY * 2 + hairpin_init[3], ComputeEnergy(parsing::ParseViennaRna("AAAAU", "(...)")));
}

}
}
