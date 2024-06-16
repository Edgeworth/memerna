// Copyright 2024 Eliot Courtney.
#include "model/energy.h"

#include <cmath>
#include <cstring>
#include <string>

#include "gtest/gtest.h"
#include "tests/init.h"
#include "tests/util.h"

namespace mrna {

class EnergyTestT04 : public testing::TestWithParam<int> {};

#if ENERGY_PRECISION == 1

TEST_P(EnergyTestT04, T04P1) {
  auto m = t04_ms[GetParam()];

  EXPECT_EQ(E(4.5), GetEnergy(m, "GCAAAGCC", "((...).)"));
  EXPECT_EQ(E(5.7), GetEnergy(m, "CCCAAAAUG", ".(.(...))"));
  EXPECT_EQ(E(5.5), GetEnergy(m, "UACAGA", "(....)"));
  EXPECT_EQ(E(-0.6), GetEnergy(m, "AGGGUCAUCCG", ".(((...)))."));
  EXPECT_EQ(E(8.0), GetEnergy(m, "AGAGAAACAAAU", "(..(...)...)"));
  EXPECT_EQ(E(9.5), GetEnergy(m, "CGUUGCCUAAAAAGGAAACAAG", "(.............(...)..)"));
  EXPECT_EQ(E(7.7), GetEnergy(m, "CCCGAAACAG", "(..(...).)"));
  EXPECT_EQ(E(7.4), GetEnergy(m, "GACAGAAACGCUGAAUC", "((..(...)......))"));
  EXPECT_EQ(E(17.3), GetEnergy(m, "CUGAAACUGGAAACAGAAAUG", "(.(...)..(...).(...))"));
  EXPECT_EQ(E(18.2), GetEnergy(m, "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))"));
  EXPECT_EQ(E(17.6), GetEnergy(m, "AGCUAAAAACAAAGGUGAAACGU", "(..(...).(...)..(...).)"));
  EXPECT_EQ(E(13.1), GetEnergy(m, "CUGAAACUGGAAACAGAAAUG", ".(.(...)(....)......)"));
  EXPECT_EQ(E(-27.6),
      GetEnergy(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
          "(((((((((((.((...((((....))))..)).)))..((((..((((....))))...)))).))))))))...."));
  EXPECT_EQ(E(17.9), GetEnergy(m, "UCUGAGUAAAUUGCUACGCG", "(....)((...).......)"));
  EXPECT_EQ(E(-43.1), GetEnergy(m, k16sHSapiens3));

  // Special stacking - this is not implemented. TODO(4): Implement this?
  EXPECT_EQ(E(3.7), GetEnergy(m, "GGUCAAAGGUC", "((((...))))"));
  EXPECT_EQ(E(-4.5), GetEnergy(m, "GGGGAAACCCC", "((((...))))"));
  EXPECT_EQ(E(7.2), GetEnergy(m, "UGACAAAGGCGA", "(..(...)...)"));
}

#elif ENERGY_PRECISION == 2

TEST_P(EnergyTestT04, T04P2) {
  auto m = t04_ms[GetParam()];

  EXPECT_EQ(E(4.45), GetEnergy(m, "GCAAAGCC", "((...).)"));
  EXPECT_EQ(E(5.71), GetEnergy(m, "CCCAAAAUG", ".(.(...))"));
  EXPECT_EQ(E(5.50), GetEnergy(m, "UACAGA", "(....)"));
  EXPECT_EQ(E(-0.59), GetEnergy(m, "AGGGUCAUCCG", ".(((...)))."));
  EXPECT_EQ(E(8.00), GetEnergy(m, "AGAGAAACAAAU", "(..(...)...)"));
  EXPECT_EQ(E(9.50), GetEnergy(m, "CGUUGCCUAAAAAGGAAACAAG", "(.............(...)..)"));
  EXPECT_EQ(E(7.70), GetEnergy(m, "CCCGAAACAG", "(..(...).)"));
  EXPECT_EQ(E(7.45), GetEnergy(m, "GACAGAAACGCUGAAUC", "((..(...)......))"));
  EXPECT_EQ(E(17.30), GetEnergy(m, "CUGAAACUGGAAACAGAAAUG", "(.(...)..(...).(...))"));
  EXPECT_EQ(E(18.25), GetEnergy(m, "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))"));
  EXPECT_EQ(E(17.60), GetEnergy(m, "AGCUAAAAACAAAGGUGAAACGU", "(..(...).(...)..(...).)"));
  EXPECT_EQ(E(13.10), GetEnergy(m, "CUGAAACUGGAAACAGAAAUG", ".(.(...)(....)......)"));
  EXPECT_EQ(E(-27.41),
      GetEnergy(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
          "(((((((((((.((...((((....))))..)).)))..((((..((((....))))...)))).))))))))...."));
  EXPECT_EQ(E(17.90), GetEnergy(m, "UCUGAGUAAAUUGCUACGCG", "(....)((...).......)"));
  EXPECT_EQ(E(-42.75), GetEnergy(m, k16sHSapiens3));

  // Special stacking - this is not implemented.
  EXPECT_EQ(E(3.63), GetEnergy(m, "GGUCAAAGGUC", "((((...))))"));
  EXPECT_EQ(E(-4.38), GetEnergy(m, "GGGGAAACCCC", "((((...))))"));
  EXPECT_EQ(E(7.20), GetEnergy(m, "UGACAAAGGCGA", "(..(...)...)"));
}

#endif

INSTANTIATE_TEST_SUITE_P(EnergyTest, EnergyTestT04, testing::Range(0, NUM_T04_MODELS));

}  // namespace mrna
