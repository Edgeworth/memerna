// Copyright 2022 Eliot Courtney.
#include "model/model.h"

#include <memory>

#include "gtest/gtest.h"

namespace mrna {

TEST(ModelTest, EnergyToString) {
  EXPECT_EQ("0.00", EnergyToString(0));
  EXPECT_EQ("0.01", EnergyToString(1));
  EXPECT_EQ("0.10", EnergyToString(10));
  EXPECT_EQ("1.00", EnergyToString(100));
  EXPECT_EQ("10.00", EnergyToString(1000));
  EXPECT_EQ("-0.01", EnergyToString(-1));
  EXPECT_EQ("-0.10", EnergyToString(-10));
  EXPECT_EQ("-1.00", EnergyToString(-100));
  EXPECT_EQ("-10.00", EnergyToString(-1000));
  EXPECT_EQ("MAX", EnergyToString(MAX_E));
  EXPECT_EQ("CAP", EnergyToString(CAP_E));
}

TEST(ModelTest, EnergyFromString) {
  EXPECT_EQ(0, EnergyFromString("0.00"));
  EXPECT_EQ(1, EnergyFromString("0.01"));
  EXPECT_EQ(10, EnergyFromString("0.10"));
  EXPECT_EQ(100, EnergyFromString("1.00"));
  EXPECT_EQ(1000, EnergyFromString("10.00"));
  EXPECT_EQ(-1, EnergyFromString("-0.01"));
  EXPECT_EQ(-10, EnergyFromString("-0.10"));
  EXPECT_EQ(-100, EnergyFromString("-1.00"));
  EXPECT_EQ(-1000, EnergyFromString("-10.00"));
  EXPECT_EQ(MAX_E, EnergyFromString("MAX"));
  EXPECT_EQ(CAP_E, EnergyFromString("CAP"));
}

}  // namespace mrna
