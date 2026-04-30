// Copyright 2026 Eliot Courtney.
#include "api/energy/energy_cfg.h"

#include <stdexcept>
#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "util/argparse.h"

namespace mrna::erg {

namespace {

std::vector<char*> MakeArgv(std::vector<std::string>& args) {
  std::vector<char*> argv;
  argv.reserve(args.size());
  for (auto& s : args) argv.push_back(s.data());
  return argv;
}

ArgParse ParseEnergyArgs(std::vector<std::string> argv_str) {
  ArgParse args;
  RegisterOptsEnergyCfg(&args);
  auto argv = MakeArgv(argv_str);
  EXPECT_EQ("", args.Parse(static_cast<int>(argv.size()), argv.data()));
  return args;
}

}  // namespace

TEST(EnergyCfgTest, PfnFromArgParseNormalizesDefaultBulgeStates) {
  const auto cfg = EnergyCfgPfn::FromArgParse(ParseEnergyArgs({"prog"}));
  EXPECT_FALSE(cfg.bulge_states);
}

TEST(EnergyCfgTest, PfnFromArgParseAllowsExplicitNoBulgeStates) {
  const auto cfg = EnergyCfgPfn::FromArgParse(ParseEnergyArgs({"prog", "--no-bulge-states"}));
  EXPECT_FALSE(cfg.bulge_states);
}

TEST(EnergyCfgTest, PfnFromArgParseRejectsExplicitBulgeStates) {
  const auto args = ParseEnergyArgs({"prog", "--bulge-states"});
  EXPECT_THROW((void)EnergyCfgPfn::FromArgParse(args), std::runtime_error);
}

}  // namespace mrna::erg
