// Copyright 2026 Eliot Courtney.
#include "util/argparse.h"

#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "tests/util.h"
#include "util/log.h"

namespace mrna {
namespace {

void ParseHelp(const char* help) {
  const auto required = mrna::Opt(mrna::Opt::ARG).LongName("required").Required();
  mrna::ArgParse args;
  args.RegisterOpt(required);
  std::vector<std::string> argv_str = {"prog", help};
  auto argv = MakeArgv(argv_str);
  args.ParseOrExit(static_cast<int>(argv.size()), argv.data());
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

TEST(ArgParseTest, StringArgumentPreservesSpaces) {
  const auto opt = mrna::Opt(mrna::Opt::ARG).LongName("path");
  mrna::ArgParse args;
  args.RegisterOpt(opt);

  std::vector<std::string> argv_str = {"prog", "--path", "/tmp/model path"};
  auto argv = MakeArgv(argv_str);
  EXPECT_EQ("", args.Parse(static_cast<int>(argv.size()), argv.data()));

  EXPECT_EQ("/tmp/model path", args.Get<std::string>(opt));
}

TEST(ArgParseTest, StringPositionalPreservesSpaces) {
  mrna::ArgParse args;

  std::vector<std::string> argv_str = {"prog", "RNA with spaces"};
  auto argv = MakeArgv(argv_str);
  EXPECT_EQ("", args.Parse(static_cast<int>(argv.size()), argv.data()));

  EXPECT_EQ("RNA with spaces", args.Pos<std::string>(0));
}

TEST(ArgParseTest, ShortOnlyFlagParses) {
  const auto flag = mrna::Opt(mrna::Opt::FLAG).ShortName("x");
  mrna::ArgParse args;
  args.RegisterOpt(flag);

  std::vector<std::string> argv_str = {"prog", "-x"};
  auto argv = MakeArgv(argv_str);
  EXPECT_EQ("", args.Parse(static_cast<int>(argv.size()), argv.data()));

  EXPECT_TRUE(args.HasExplicit(flag));
  EXPECT_TRUE(args.Get<bool>(flag));
}

TEST(ArgParseTest, RegisteringTheSameOptionTwiceDoesNotDuplicateUsage) {
  const auto opt = mrna::Opt(mrna::Opt::FLAG).LongName("flag");
  mrna::ArgParse args;
  args.RegisterOpt(opt);
  args.RegisterOpt(opt);

  const auto usage = args.Usage();
  const auto first = usage.find("--flag");
  ASSERT_NE(std::string::npos, first);
  EXPECT_EQ(std::string::npos, usage.find("--flag", first + 1));
}

TEST(ArgParseTest, UsageIncludesProgramAndHelp) {
  mrna::ArgParse args;

  const auto usage = args.Usage("prog");

  EXPECT_NE(std::string::npos, usage.find("Usage: prog"));
  EXPECT_NE(std::string::npos, usage.find("-h, --help"));
}

TEST(ArgParseTest, ParseOrExitHelpSucceedsBeforeRequiredValidation) {
  EXPECT_EXIT(ParseHelp("--help"), testing::ExitedWithCode(0), "");
  EXPECT_EXIT(ParseHelp("-h"), testing::ExitedWithCode(0), "");
}

}  // namespace mrna
