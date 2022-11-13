// Copyright 2016 Eliot Courtney.
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <deque>
#include <random>
#include <string>
#include <utility>

#include "bridge/bridge.h"
#include "fuzz/fuzz_cfg.h"
#include "fuzz/fuzz_invocation.h"
#include "model/primary.h"
#include "programs/fuzz/fuzz_harness.h"
#include "util/argparse.h"
#include "util/error.h"

inline const auto OPT_PRINT_INTERVAL = mrna::Opt(mrna::Opt::ARG)
                                           .LongName("print-interval")
                                           .Default("5")
                                           .Help("status update every n seconds");
inline const auto OPT_ENUMERATE =
    mrna::Opt(mrna::Opt::FLAG).LongName("enumerate").Help("enumerate all sequences");

int main(int argc, char* argv[]) {
  mrna::ArgParse args;
  mrna::fuzz::RegisterOpts(&args);
  args.RegisterOpt(mrna::bridge::OPT_RNASTRUCTURE_DATA);
  args.RegisterOpt(OPT_PRINT_INTERVAL);
  args.RegisterOpt(OPT_RANDOM_MODEL);
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
    error("require min and max length or a sequence");
  }

  auto harness = FuzzHarness(std::move(args));
  if (!seq.empty()) {
    const auto r = mrna::Primary::FromSeq(seq);
    printf("Running single fuzz on %s\n", seq.c_str());
    auto invoc = harness.CreateInvocation(r);
    const auto res = invoc.Run();
    if (!res.empty()) {
      for (const auto& s : res) printf("%s\n", s.c_str());
      printf("\n");
    }
  } else if (enumerate) {
    printf("Exhaustive fuzzing [%d, %d] len RNAs\n", min_len, max_len);
    int64_t i = 0;
    mrna::Primary r(min_len);
    while (static_cast<int>(r.size()) <= max_len) {
      if (++i % 1000 == 0)
        printf("Fuzzed %" PRId64 " RNA, current size: %d\n", i, static_cast<int>(r.size()));

      auto invoc = harness.CreateInvocation(r);
      const auto res = invoc.Run();
      if (!res.empty()) {
        for (const auto& s : res) printf("%s\n", s.c_str());
        printf("\n");
      }

      r.Increment();
    }
    printf("Finished exhaustive fuzzing [%d, %d] len RNAs\n", min_len, max_len);
  } else {
    printf("Random fuzzing [%d, %d] len RNAs\n", min_len, max_len);
    std::uniform_int_distribution<int> len_dist(min_len, max_len);
    auto start_time = std::chrono::steady_clock::now();
    for (int64_t i = 0;; ++i) {
      if (interval > 0 &&
          std::chrono::duration_cast<std::chrono::seconds>(
              std::chrono::steady_clock::now() - start_time)
                  .count() > interval) {
        printf("Fuzzed %" PRId64 " RNA\n", i);
        start_time = std::chrono::steady_clock::now();
      }
      int len = len_dist(harness.e());
      auto r = mrna::Primary::Random(len);

      auto invoc = harness.CreateInvocation(r);
      const auto res = invoc.Run();
      if (!res.empty()) {
        for (const auto& s : res) printf("%s\n", s.c_str());
        printf("\n");
      }
    }
  }
}
