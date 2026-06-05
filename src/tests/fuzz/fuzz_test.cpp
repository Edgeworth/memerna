// Copyright 2026 Eliot Courtney.
#include <stdexcept>

#include "fuzz/fuzz_cfg.h"
#include "gtest/gtest.h"
#include "tests/util.h"

namespace mrna::fuzz {

TEST(FuzzCfgTest, PfnSuboptAllowsDefaultBulgeStates) {
  const auto cfg =
      FuzzCfg::FromArgParse(ParseTestArgs({"prog", "--pfn-subopt", "8"}, RegisterOpts));
  EXPECT_EQ(8, cfg.pfn_subopt);
  EXPECT_TRUE(cfg.pfn);
  EXPECT_TRUE(cfg.energy_cfg.bulge_states);
}

TEST(FuzzCfgTest, PfnSuboptAllowsExplicitNoBulgeStates) {
  const auto cfg = FuzzCfg::FromArgParse(
      ParseTestArgs({"prog", "--pfn-subopt", "8", "--no-bulge-states"}, RegisterOpts));
  EXPECT_EQ(8, cfg.pfn_subopt);
  EXPECT_TRUE(cfg.pfn);
  EXPECT_FALSE(cfg.energy_cfg.bulge_states);
}

TEST(FuzzCfgTest, PfnSuboptRejectsExplicitBulgeStates) {
  const auto args = ParseTestArgs({"prog", "--pfn-subopt", "8", "--bulge-states"}, RegisterOpts);
  EXPECT_THROW((void)FuzzCfg::FromArgParse(args), std::runtime_error);
}

}  // namespace mrna::fuzz
