// Copyright 2026 Eliot Courtney.
#include "util/argparse.h"

#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "util/log.h"

namespace mrna {

namespace {

std::vector<char*> MakeArgv(std::vector<std::string>& args) {
  std::vector<char*> argv;
  argv.reserve(args.size());
  for (auto& s : args) argv.push_back(s.data());
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
  const auto old_level = GetLogLevel();
  SetLogLevel(LogLevel::INFO);

  mrna::ArgParse args;
  std::vector<std::string> argv_str = {"prog", "--no-verbose"};
  auto argv = MakeArgv(argv_str);
  args.ParseOrExit(static_cast<int>(argv.size()), argv.data());

  EXPECT_NE(GetLogLevel(), LogLevel::DEBUG);
  SetLogLevel(old_level);
}

TEST(ArgParseTest, ParseOrExitVerboseEnablesDebug) {
  const auto old_level = GetLogLevel();
  SetLogLevel(LogLevel::INFO);

  mrna::ArgParse args;
  std::vector<std::string> argv_str = {"prog", "--verbose"};
  auto argv = MakeArgv(argv_str);
  args.ParseOrExit(static_cast<int>(argv.size()), argv.data());

  EXPECT_EQ(GetLogLevel(), LogLevel::DEBUG);
  SetLogLevel(old_level);
}

TEST(ArgParseTest, ThreeDashesRejected) {
  const auto opt = mrna::Opt(mrna::Opt::FLAG).LongName("verbose");
  mrna::ArgParse args;
  args.RegisterOpt(opt);

  std::vector<std::string> argv_str = {"prog", "---verbose"};
  auto argv = MakeArgv(argv_str);
  EXPECT_NE("", args.Parse(static_cast<int>(argv.size()), argv.data()));
}

TEST(ArgParseTest, BareDashIsPositional) {
  mrna::ArgParse args;

  std::vector<std::string> argv_str = {"prog", "-"};
  auto argv = MakeArgv(argv_str);
  EXPECT_EQ("", args.Parse(static_cast<int>(argv.size()), argv.data()));
  EXPECT_EQ(1, args.Pos().size());
  EXPECT_EQ("-", args.Pos()[0]);
}

}  // namespace mrna
