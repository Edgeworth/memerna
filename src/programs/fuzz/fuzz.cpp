// Copyright 2016 Eliot Courtney.
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <deque>
#include <random>
#include <string>
#include <utility>

#include "bridge/bridge.h"
#include "fuzz/config.h"
#include "fuzz/fuzz_invocation.h"
#include "model/primary.h"
#include "programs/fuzz/fuzz_harness.h"
#include "util/argparse.h"
#include "util/error.h"

inline const auto OPT_PRINT_INTERVAL = mrna::Opt(mrna::Opt::ARG)
                                           .LongName("print-interval")
                                           .Default("5")
                                           .Help("status update every n seconds");

int main(int argc, char* argv[]) {
  mrna::ArgParse args;
  mrna::fuzz::RegisterOpts(&args);
  args.RegisterOpt(mrna::bridge::OPT_RNASTRUCTURE_DATA);
  args.RegisterOpt(OPT_PRINT_INTERVAL);
  args.RegisterOpt(OPT_RANDOM);
  args.ParseOrExit(argc, argv);

  verify(args.PosSize() == 2, "require min and max length");
  const int min_len = args.Pos<int>(0);
  const int max_len = args.Pos<int>(1);
  const auto interval = args.Get<int>(OPT_PRINT_INTERVAL);
  verify(min_len > 0, "invalid min length");
  verify(max_len >= min_len, "invalid max len");
  std::uniform_int_distribution<int> len_dist(min_len, max_len);

  auto harness = FuzzHarness(std::move(args));

  printf("Regular fuzzing [%d, %d] len RNAs\n", min_len, max_len);
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
