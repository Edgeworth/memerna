// Copyright 2022 Eliot Courtney.
#include "gtest/gtest.h"
#include "model/energy.h"

namespace mrna {

#if ENERGY_PRECISION == 1

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

#elif ENERGY_PRECISION == 2

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

}  // namespace mrna
