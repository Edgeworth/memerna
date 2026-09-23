// Copyright 2026 Eliot Courtney.
#include <stdexcept>
#include <string>
#include <vector>

#include "fuzz/fuzz_cfg.h"
#include "gtest/gtest.h"
#include "tests/util.h"

namespace mrna::fuzz {

TEST(FuzzCfgTest, PfnSuboptAllowsDefaultBulgeStates) {
  const auto cfg =
      FuzzCfg::FromArgParse(ParseTestArgs({"prog", "--pfn-subopt", "8"}, RegisterOpts));
  EXPECT_EQ(8, cfg.pfn_subopt);
  EXPECT_FALSE(cfg.pfn);
  EXPECT_TRUE(cfg.PfnEnabled());
  EXPECT_TRUE(cfg.energy_cfg.bulge_states);
}

TEST(FuzzCfgTest, PfnSuboptAllowsExplicitNoBulgeStates) {
  const auto cfg = FuzzCfg::FromArgParse(
      ParseTestArgs({"prog", "--pfn-subopt", "8", "--no-bulge-states"}, RegisterOpts));
  EXPECT_EQ(8, cfg.pfn_subopt);
  EXPECT_FALSE(cfg.pfn);
  EXPECT_TRUE(cfg.PfnEnabled());
  EXPECT_FALSE(cfg.energy_cfg.bulge_states);
}

TEST(FuzzCfgTest, PfnSuboptRejectsExplicitBulgeStates) {
  const auto args = ParseTestArgs({"prog", "--pfn-subopt", "8", "--bulge-states"}, RegisterOpts);
  EXPECT_THROW((void)FuzzCfg::FromArgParse(args), std::runtime_error);
}

#ifdef MRNA_USE_RNASTRUCTURE
TEST(FuzzCfgTest, RejectsUnsupportedRNAstructureModel) {
  const auto args = ParseTestArgs(
      {"prog", "--mfe-rnastructure", "--energy-model", "t12", "--rnastructure-data", "data"},
      RegisterOpts);
  EXPECT_THROW((void)FuzzCfg::FromArgParse(args), std::runtime_error);
}

TEST(FuzzCfgTest, RejectsUnsupportedRNAstructureConfigurations) {
  const std::vector<std::vector<std::string>> cases{
      {"prog", "--mfe-rnastructure", "--random-seeds"},
      {"prog", "--subopt-rnastructure", "--random-pf"},
      {"prog", "--pfn-rnastructure", "--lonely-pairs", "off"},
      {"prog", "--pfn-rnastructure", "--ctd", "none"},
      {"prog", "--mfe-rnastructure", "--no-bulge-states"},
      {"prog", "--subopt-rnastructure", "--no-bulge-states"},
  };
  for (auto argv : cases) {
    SCOPED_TRACE(argv[1]);
    argv.insert(argv.end(), {"--rnastructure-data", "data"});
    const auto args = ParseTestArgs(argv, RegisterOpts);
    EXPECT_THROW((void)FuzzCfg::FromArgParse(args), std::runtime_error);
  }
}

TEST(FuzzCfgTest, RNAstructurePfnAllowsNoBulgeStates) {
  const auto cfg = FuzzCfg::FromArgParse(ParseTestArgs(
      {"prog", "--pfn-rnastructure", "--no-bulge-states", "--rnastructure-data", "data"},
      RegisterOpts));
  EXPECT_TRUE(cfg.pfn_rnastructure);
  EXPECT_FALSE(cfg.energy_cfg.bulge_states);
}
#else
TEST(FuzzCfgTest, RejectsRNAstructureChecksWhenUnavailable) {
  const auto args = ParseTestArgs({"prog", "--mfe-rnastructure"}, RegisterOpts);
  EXPECT_THROW((void)FuzzCfg::FromArgParse(args), std::runtime_error);
}
#endif  // MRNA_USE_RNASTRUCTURE

}  // namespace mrna::fuzz
