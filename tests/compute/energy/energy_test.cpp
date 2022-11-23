// Copyright 2022 Eliot Courtney.
#include "compute/energy/energy.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "common_test.h"
#include "compute/energy/branch.h"
#include "compute/energy/t04/precomp.h"
#include "gtest/gtest.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/primary.h"
#include "model/secondary.h"

namespace mrna::energy {

TEST(EnergyTest, Helpers) {
  EXPECT_EQ(0, MaxNumContiguous(Primary::FromSeq("")));
  EXPECT_EQ(1, MaxNumContiguous(Primary::FromSeq("A")));
  EXPECT_EQ(2, MaxNumContiguous(Primary::FromSeq("AA")));
  EXPECT_EQ(2, MaxNumContiguous(Primary::FromSeq("GUAAC")));
  EXPECT_EQ(1, MaxNumContiguous(Primary::FromSeq("GUACA")));
  EXPECT_EQ(3, MaxNumContiguous(Primary::FromSeq("GAUCCC")));
  EXPECT_EQ(3, MaxNumContiguous(Primary::FromSeq("GGGAUC")));
  EXPECT_EQ(4, MaxNumContiguous(Primary::FromSeq("GGGAUCAAAA")));
  EXPECT_EQ(5, MaxNumContiguous(Primary::FromSeq("GGGAUUUUUCAAAA")));
}

TEST(EnergyTest, GetBranchCounts) {
  EXPECT_EQ((std::vector<int>{2, 0}), GetBranchCounts(Secondary::FromDb("()")));
  EXPECT_EQ((std::vector<int>{}), GetBranchCounts(Secondary::FromDb("")));
  EXPECT_EQ((std::vector<int>{0}), GetBranchCounts(Secondary::FromDb(".")));
  EXPECT_EQ((std::vector<int>{2, 0, 2, 0, 2, 0}), GetBranchCounts(Secondary::FromDb("()()()")));
  EXPECT_EQ((std::vector<int>{2, 1, 0, 1}), GetBranchCounts(Secondary::FromDb("(())")));
  EXPECT_EQ(
      (std::vector<int>{2, 3, 0, 3, 0, 3, 0, 3}), GetBranchCounts(Secondary::FromDb("(()()())")));
  EXPECT_EQ(
      (std::vector<int>{2, 1, 2, 0, 2, 0, 2, 1}), GetBranchCounts(Secondary::FromDb("((()()))")));
}

}  // namespace mrna::energy
