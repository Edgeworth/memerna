// Copyright 2022 Eliot Courtney.
#include <stdexcept>

#include "api/energy/pseudofree_cfg.h"
#include "backends/common/energy.h"
#include "gtest/gtest.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"

namespace mrna {

#if MRNA_ENERGY_PRECISION == 1

TEST(ModelTest, EnergyToString) {
  EXPECT_EQ("0.0", E(0).ToString());
  EXPECT_EQ("0.1", E(0.1).ToString());
  EXPECT_EQ("1.0", E(1.0).ToString());
  EXPECT_EQ("10.0", E(10.0).ToString());
  EXPECT_EQ("-0.1", E(-0.1).ToString());
  EXPECT_EQ("-1.0", E(-1.0).ToString());
  EXPECT_EQ("-10.0", E(-10.0).ToString());
  EXPECT_EQ("MAX", MAX_E.ToString());
  EXPECT_EQ("CAP", CAP_E.ToString());
}

TEST(ModelTest, EnergyFromString) {
  EXPECT_EQ(E(0.0), Energy::FromString("0.0"));
  EXPECT_EQ(E(0.1), Energy::FromString("0.1"));
  EXPECT_EQ(E(1.0), Energy::FromString("1.0"));
  EXPECT_EQ(E(10.0), Energy::FromString("10.0"));
  EXPECT_EQ(E(-0.1), Energy::FromString("-0.1"));
  EXPECT_EQ(E(-1.0), Energy::FromString("-1.0"));
  EXPECT_EQ(E(-10.0), Energy::FromString("-10.0"));
  EXPECT_EQ(MAX_E, Energy::FromString("MAX"));
  EXPECT_EQ(CAP_E, Energy::FromString("CAP"));
}

TEST(ModelTest, EnergyFromRaw) {
  EXPECT_EQ(Energy::FromRaw(0), E(0.00));
  EXPECT_EQ(Energy::FromRaw(1), E(0.1));
  EXPECT_EQ(Energy::FromRaw(10), E(1.0));
  EXPECT_EQ(Energy::FromRaw(100), E(10.0));
  EXPECT_EQ(Energy::FromRaw(1000), E(100.0));
  EXPECT_EQ(Energy::FromRaw(-1), E(-0.1));
  EXPECT_EQ(Energy::FromRaw(-10), E(-1.0));
  EXPECT_EQ(Energy::FromRaw(-100), E(-10.0));
  EXPECT_EQ(Energy::FromRaw(-1000), E(-100.0));
}

#elif MRNA_ENERGY_PRECISION == 2

TEST(ModelTest, EnergyToString) {
  EXPECT_EQ("0.00", E(0).ToString());
  EXPECT_EQ("0.01", E(0.01).ToString());
  EXPECT_EQ("0.10", E(0.10).ToString());
  EXPECT_EQ("1.00", E(1.00).ToString());
  EXPECT_EQ("10.00", E(10.00).ToString());
  EXPECT_EQ("-0.01", E(-0.01).ToString());
  EXPECT_EQ("-0.10", E(-0.10).ToString());
  EXPECT_EQ("-1.00", E(-1.00).ToString());
  EXPECT_EQ("-10.00", E(-10.00).ToString());
  EXPECT_EQ("MAX", MAX_E.ToString());
  EXPECT_EQ("CAP", CAP_E.ToString());
}

TEST(ModelTest, EnergyFromString) {
  EXPECT_EQ(E(0.00), Energy::FromString("0.00"));
  EXPECT_EQ(E(0.01), Energy::FromString("0.01"));
  EXPECT_EQ(E(0.10), Energy::FromString("0.10"));
  EXPECT_EQ(E(1.00), Energy::FromString("1.00"));
  EXPECT_EQ(E(10.00), Energy::FromString("10.00"));
  EXPECT_EQ(E(-0.01), Energy::FromString("-0.01"));
  EXPECT_EQ(E(-0.10), Energy::FromString("-0.10"));
  EXPECT_EQ(E(-1.00), Energy::FromString("-1.00"));
  EXPECT_EQ(E(-10.00), Energy::FromString("-10.00"));
  EXPECT_EQ(MAX_E, Energy::FromString("MAX"));
  EXPECT_EQ(CAP_E, Energy::FromString("CAP"));
}

TEST(ModelTest, EnergyFromRaw) {
  EXPECT_EQ(Energy::FromRaw(0), E(0.00));
  EXPECT_EQ(Energy::FromRaw(1), E(0.01));
  EXPECT_EQ(Energy::FromRaw(10), E(0.10));
  EXPECT_EQ(Energy::FromRaw(100), E(1.00));
  EXPECT_EQ(Energy::FromRaw(1000), E(10.00));
  EXPECT_EQ(Energy::FromRaw(-1), E(-0.01));
  EXPECT_EQ(Energy::FromRaw(-10), E(-0.10));
  EXPECT_EQ(Energy::FromRaw(-100), E(-1.00));
  EXPECT_EQ(Energy::FromRaw(-1000), E(-10.00));
}

#endif

TEST(ModelTest, MinEnergyTraversesNestedArrays) {
  const Energy energies[2][2] = {{E(4.0), E(3.0)}, {E(2.0), E(-1.0)}};

  EXPECT_EQ(E(-1.0), md::MinEnergy(energies));
}

TEST(PseudofreeCfgTest, BoltzRangesHandleCappedEntries) {
  const auto pf = erg::PseudofreeCfg({}, {CAP_E, ZERO_E}).Boltz();

  EXPECT_EQ(pf.UnpairedProd(0, 0), ZERO_B);
  EXPECT_EQ(pf.UnpairedProd(0, 1), ZERO_B);
  EXPECT_EQ(pf.UnpairedProd(1, 1), ONE_B);
  EXPECT_EQ(pf.UnpairedProd(1, 0), ONE_B);
}

TEST(SecondaryTest, VerifyChecksPairsAgainstSequence) {
  EXPECT_NO_THROW(Secondary::FromDb("((...)).").Verify(Primary::FromSeq("GGAAACCA")));
  // G-U is a canonical pair.
  EXPECT_NO_THROW(Secondary::FromDb("(...)").Verify(Primary::FromSeq("GAAAU")));
  // G-A is not.
  EXPECT_THROW(Secondary::FromDb("(...)").Verify(Primary::FromSeq("GAAAA")), std::runtime_error);
  // The hairpin has 2 unpaired bases.
  EXPECT_THROW(Secondary::FromDb("(..)").Verify(Primary::FromSeq("GAAC")), std::runtime_error);
  // The lengths differ.
  EXPECT_THROW(Secondary::FromDb("(...)").Verify(Primary::FromSeq("GAAACA")), std::runtime_error);
}

}  // namespace mrna
