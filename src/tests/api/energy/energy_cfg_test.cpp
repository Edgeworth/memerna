// Copyright 2026 Eliot Courtney.
#include "api/energy/energy_cfg.h"

#include <stdexcept>

#include "gtest/gtest.h"
#include "tests/util.h"

namespace mrna::erg {

TEST(EnergyCfgTest, PfnFromArgParseNormalizesDefaultBulgeStates) {
  const auto cfg = EnergyCfgPfn::FromArgParse(ParseTestArgs({"prog"}, RegisterOptsEnergyCfg));
  EXPECT_FALSE(cfg.bulge_states);
}

TEST(EnergyCfgTest, PfnFromArgParseAllowsExplicitNoBulgeStates) {
  const auto cfg = EnergyCfgPfn::FromArgParse(
      ParseTestArgs({"prog", "--no-bulge-states"}, RegisterOptsEnergyCfg));
  EXPECT_FALSE(cfg.bulge_states);
}

TEST(EnergyCfgTest, PfnFromArgParseRejectsExplicitBulgeStates) {
  const auto args = ParseTestArgs({"prog", "--bulge-states"}, RegisterOptsEnergyCfg);
  EXPECT_THROW((void)EnergyCfgPfn::FromArgParse(args), std::runtime_error);
}

}  // namespace mrna::erg
