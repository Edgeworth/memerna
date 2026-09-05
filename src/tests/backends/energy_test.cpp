// Copyright 2022 Eliot Courtney.

#include <stdexcept>
#include <variant>
#include <vector>

#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "gtest/gtest.h"
#include "model/branch.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "tests/init.h"

namespace mrna::md {

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

TEST(EnergyTest, RejectsEmptySequences) {
  const Primary r;
  const Secondary s;
  for (const auto& model : t04_ms) {
    std::visit(
        [&](const auto& m) {
          EXPECT_THROW(
              (void)m->template TotalEnergy<erg::EnergyCfg{}>(r, s, nullptr, erg::PseudofreeCfg{}),
              std::runtime_error);
        },
        model);
  }
}

}  // namespace mrna::md
