// Copyright 2016 Eliot Courtney.
#include <fmt/core.h>

#include <chrono>
#include <cinttypes>
#include <deque>
#include <ios>
#include <random>
#include <string>
#include <utility>

#include "api/bridge/bridge.h"
#include "fuzz/fuzz_cfg.h"
#include "fuzz/fuzz_harness.h"
#include "fuzz/fuzz_invocation.h"
#include "model/primary.h"
#include "util/argparse.h"
#include "util/error.h"

inline const auto OPT_PRINT_INTERVAL = mrna::Opt(mrna::Opt::ARG)
                                           .LongName("print-interval")
                                           .Default("5")
                                           .Help("status update every n seconds");
inline const auto OPT_ENUMERATE =
    mrna::Opt(mrna::Opt::FLAG).LongName("enumerate").Help("enumerate all sequences");

int main(int argc, char* argv[]) {
  mrna::InitProgram();
  mrna::ArgParse args;
  mrna::fuzz::RegisterOpts(&args);
  args.RegisterOpt(OPT_PRINT_INTERVAL);
  args.RegisterOpt(OPT_ENUMERATE);
  args.ParseOrExit(argc, argv);

  const auto interval = args.Get<int>(OPT_PRINT_INTERVAL);
  const auto enumerate = args.Has(OPT_ENUMERATE);
  int min_len = 0;
  int max_len = 0;
  std::string seq;

  if (args.PosSize() == 1) {
    seq = args.Pos(0);
  } else if (args.PosSize() == 2) {
    min_len = args.Pos<int>(0);
    max_len = args.Pos<int>(1);
    verify(min_len > 0, "invalid min length");
    verify(max_len >= min_len, "invalid max len");
  } else {
    fatal("require min and max length or a sequence");
  }

  auto fuzz_cfg = mrna::fuzz::FuzzCfg::FromArgParse(args);
  auto harness = mrna::fuzz::FuzzHarness(fuzz_cfg);
  if (!seq.empty()) {
    const auto r = mrna::Primary::FromSeq(seq);
    fmt::print("Running single fuzz on {}\n", seq);
    auto invoc = harness.CreateInvocation(r);
    const auto res = invoc.Run();
    if (!res.empty()) {
      for (const auto& s : res) fmt::print("{}\n", s);
      fmt::print("\n");
    }
  } else if (enumerate) {
    fmt::print("Exhaustive fuzzing [{}-{}] len RNAs\n", min_len, max_len);
    int64_t i = 0;
    mrna::Primary r(min_len);
    while (static_cast<int>(r.size()) <= max_len) {
      if (++i % 1000 == 0) fmt::print("Fuzzed {} RNA, current size: {}\n", i, r.size());

      auto invoc = harness.CreateInvocation(r);
      const auto res = invoc.Run();
      if (!res.empty()) {
        for (const auto& s : res) fmt::print("{}\n", s);
        fmt::print("\n");
      }

      r.Increment();
    }
    fmt::print("Finished exhaustive fuzzing [{}-{}] len RNAs\n", min_len, max_len);
  } else {
    fmt::print("Random fuzzing [{}-{}] len RNAs\n", min_len, max_len);
    std::uniform_int_distribution<int> len_dist(min_len, max_len);
    auto start_time = std::chrono::steady_clock::now();
    for (int64_t i = 0;; ++i) {
      if (interval > 0 &&
          std::chrono::duration_cast<std::chrono::seconds>(
              std::chrono::steady_clock::now() - start_time)
                  .count() > interval) {
        fmt::print("Fuzzed {} RNAs\n", i);
        start_time = std::chrono::steady_clock::now();
      }
      const int len = len_dist(harness.e());
      auto r = mrna::Primary::Random(len);

      auto invoc = harness.CreateInvocation(r);
      const auto res = invoc.Run();
      if (!res.empty()) {
        if (fuzz_cfg.random_models) fmt::print("Random model seed: {}\n", invoc.seed());
        for (const auto& s : res) fmt::print("{}\n", s);
        fmt::print("\n");
      }
    }
  }
}
