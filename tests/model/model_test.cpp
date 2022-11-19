// Copyright 2022 Eliot Courtney.
#include <memory>

#include "gtest/gtest.h"
#include "model/constants.h"

namespace mrna {

TEST(ModelTest, EnergyToString) {
  EXPECT_EQ("0.00", E(0).ToString());
  EXPECT_EQ("0.01", E(1).ToString());
  EXPECT_EQ("0.10", E(10).ToString());
  EXPECT_EQ("1.00", E(100).ToString());
  EXPECT_EQ("10.00", E(1000).ToString());
  EXPECT_EQ("-0.01", E(-1).ToString());
  EXPECT_EQ("-0.10", E(-10).ToString());
  EXPECT_EQ("-1.00", E(-100).ToString());
  EXPECT_EQ("-10.00", E(-1000).ToString());
  EXPECT_EQ("MAX", MAX_E.ToString());
  EXPECT_EQ("CAP", CAP_E.ToString());
}

TEST(ModelTest, EnergyFromString) {
  EXPECT_EQ(Energy::FromRaw(0), Energy::FromString("0.00"));
  EXPECT_EQ(Energy::FromRaw(1), Energy::FromString("0.01"));
  EXPECT_EQ(Energy::FromRaw(10), Energy::FromString("0.10"));
  EXPECT_EQ(Energy::FromRaw(100), Energy::FromString("1.00"));
  EXPECT_EQ(Energy::FromRaw(1000), Energy::FromString("10.00"));
  EXPECT_EQ(Energy::FromRaw(-1), Energy::FromString("-0.01"));
  EXPECT_EQ(Energy::FromRaw(-10), Energy::FromString("-0.10"));
  EXPECT_EQ(Energy::FromRaw(-100), Energy::FromString("-1.00"));
  EXPECT_EQ(Energy::FromRaw(-1000), Energy::FromString("-10.00"));
  EXPECT_EQ(MAX_E, Energy::FromString("MAX"));
  EXPECT_EQ(CAP_E, Energy::FromString("CAP"));
}

}  // namespace mrna
