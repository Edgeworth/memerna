// Copyright 2024 Eliot Courtney.
#include "model/energy.h"

#include <string>

#include "gtest/gtest.h"
#include "tests/init.h"
#include "tests/util.h"

namespace mrna {

class EnergyTestT12 : public testing::TestWithParam<int> {};

#if ENERGY_PRECISION == 1

GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(EnergyTestT12);

#elif ENERGY_PRECISION == 2

TEST_P(EnergyTestT12, T12P2) {
  const auto& m = t12_ms[GetParam()];

  // Example from https://doi.org/10.1093/nar/gkac261
  EXPECT_EQ(E(-1.74), GetEnergy(m, "UGUCGAUACCCUGUCGAUA", "((((((((...))))))))"));
  EXPECT_EQ(E(-3.08), GetEnergy(m, "UAGGUCAGCCCCUGGUCUA", "((((((((...))))))))"));

  EXPECT_EQ(E(4.45), GetEnergy(m, "GCAAAGCC", "((...).)"));
  EXPECT_EQ(E(5.71), GetEnergy(m, "CCCAAAAUG", ".(.(...))"));
  EXPECT_EQ(E(5.50), GetEnergy(m, "UACAGA", "(....)"));
  EXPECT_EQ(E(-1.36), GetEnergy(m, "AGGGUCAUCCG", ".(((...)))."));
  EXPECT_EQ(E(8.00), GetEnergy(m, "AGAGAAACAAAU", "(..(...)...)"));
  EXPECT_EQ(E(9.50), GetEnergy(m, "CGUUGCCUAAAAAGGAAACAAG", "(.............(...)..)"));
  EXPECT_EQ(E(7.70), GetEnergy(m, "CCCGAAACAG", "(..(...).)"));
  EXPECT_EQ(E(7.45), GetEnergy(m, "GACAGAAACGCUGAAUC", "((..(...)......))"));
  EXPECT_EQ(E(16.30), GetEnergy(m, "CUGAAACUGGAAACAGAAAUG", "(.(...)..(...).(...))"));
  EXPECT_EQ(E(18.25), GetEnergy(m, "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))"));
  EXPECT_EQ(E(17.60), GetEnergy(m, "AGCUAAAAACAAAGGUGAAACGU", "(..(...).(...)..(...).)"));
  EXPECT_EQ(E(12.10), GetEnergy(m, "CUGAAACUGGAAACAGAAAUG", ".(.(...)(....)......)"));
  EXPECT_EQ(E(-28.29),
      GetEnergy(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
          "(((((((((((.((...((((....))))..)).)))..((((..((((....))))...)))).))))))))...."));
  EXPECT_EQ(E(15.90), GetEnergy(m, "UCUGAGUAAAUUGCUACGCG", "(....)((...).......)"));
  EXPECT_EQ(E(-41.46), GetEnergy(m, k16sHSapiens3));

  // Special stacking - this is not implemented. TODO(4): Implement this?
  EXPECT_EQ(E(2.52), GetEnergy(m, "GGUCAAAGGUC", "((((...))))"));
  EXPECT_EQ(E(-4.38), GetEnergy(m, "GGGGAAACCCC", "((((...))))"));
  EXPECT_EQ(E(7.20), GetEnergy(m, "UGACAAAGGCGA", "(..(...)...)"));
}

#endif

INSTANTIATE_TEST_SUITE_P(EnergyTest, EnergyTestT12, testing::Range(0, NUM_T12_MODELS));

}  // namespace mrna
