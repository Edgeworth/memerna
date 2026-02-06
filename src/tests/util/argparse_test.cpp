// Copyright 2026 Eliot Courtney.
#include "util/argparse.h"

#include <spdlog/spdlog.h>

#include <string>
#include <vector>

#include "gtest/gtest.h"

namespace mrna {

namespace {

std::vector<char*> MakeArgv(const std::vector<std::string>& args) {
  std::vector<char*> argv;
  argv.reserve(args.size());
  for (const auto& s : args) argv.push_back(const_cast<char*>(s.c_str()));
  return argv;
}

}  // namespace

TEST(ArgParseTest, HasIsTrueForNoFlagButGetOrIsFalse) {
  const auto opt = mrna::Opt(mrna::Opt::FLAG).LongName("flag");
  mrna::ArgParse args;
  args.RegisterOpt(opt);

  std::vector<std::string> argv_str = {"prog", "--no-flag"};
  auto argv = MakeArgv(argv_str);
  EXPECT_EQ("", args.Parse(static_cast<int>(argv.size()), argv.data()));

  EXPECT_TRUE(args.Has(opt));
  EXPECT_FALSE(args.GetOr(opt));
}

TEST(ArgParseTest, GetOrDefaultIsFalseWhenUnspecified) {
  const auto opt = mrna::Opt(mrna::Opt::FLAG).LongName("flag");
  mrna::ArgParse args;
  args.RegisterOpt(opt);

  std::vector<std::string> argv_str = {"prog"};
  auto argv = MakeArgv(argv_str);
  EXPECT_EQ("", args.Parse(static_cast<int>(argv.size()), argv.data()));

  EXPECT_FALSE(args.Has(opt));
  EXPECT_FALSE(args.GetOr(opt));
}

TEST(ArgParseTest, ParseOrExitNoVerboseDoesNotEnableDebug) {
  const auto old_level = spdlog::get_level();
  spdlog::set_level(spdlog::level::info);

  mrna::ArgParse args;
  std::vector<std::string> argv_str = {"prog", "--no-verbose"};
  auto argv = MakeArgv(argv_str);
  args.ParseOrExit(static_cast<int>(argv.size()), argv.data());

  EXPECT_NE(spdlog::get_level(), spdlog::level::debug);
  spdlog::set_level(old_level);
}

TEST(ArgParseTest, ParseOrExitVerboseEnablesDebug) {
  const auto old_level = spdlog::get_level();
  spdlog::set_level(spdlog::level::info);

  mrna::ArgParse args;
  std::vector<std::string> argv_str = {"prog", "--verbose"};
  auto argv = MakeArgv(argv_str);
  args.ParseOrExit(static_cast<int>(argv.size()), argv.data());

  EXPECT_EQ(spdlog::get_level(), spdlog::level::debug);
  spdlog::set_level(old_level);
}

}  // namespace mrna
