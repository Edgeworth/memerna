// Copyright 2016 E.
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <ctime>
#include <deque>
#include <random>
#include <string>
#include <utility>

#include "bridge/bridge.h"
#include "compute/energy/config.h"
#include "compute/energy/model.h"
#include "fuzz/config.h"
#include "fuzz/fuzzer.h"
#include "model/primary.h"
#include "options.h"
#include "util/argparse.h"
#include "util/error.h"

#ifdef USE_RNASTRUCTURE
#include "bridge/rnastructure.h"
#endif  // USE_RNASTRUCTURE

#ifdef __AFL_FUZZ_TESTCASE_LEN
#include <unistd.h>  // For __AFL_FUZZ_TESTCASE_LEN

__AFL_FUZZ_INIT();
#endif

inline const auto OPT_PRINT_INTERVAL =
    mrna::Opt().LongName("print-interval").Default("-1").Help("status update every n seconds");

// Energy model options:
// TODO: Get this to work for RNAstructure.
inline const auto OPT_RANDOM = mrna::Opt().LongName("random").Help("use random energy models");

int main(int argc, char* argv[]) {
  std::mt19937 eng(uint_fast32_t(time(nullptr)));
  mrna::ArgParse args;
  mrna::fuzz::RegisterOpts(&args);
  args.RegisterOpt(mrna::OPT_AFL);
  args.RegisterOpt(OPT_PRINT_INTERVAL);
  args.RegisterOpt(OPT_RANDOM);
  args.RegisterOpt(mrna::energy::OPT_MEMERNA_DATA);
  args.RegisterOpt(mrna::bridge::OPT_RNASTRUCTURE_DATA);

  args.ParseOrExit(argc, argv);

  // TODO: Actually set this.
#ifdef USE_RNASTRUCTURE
  mrna::bridge::RNAstructure rnastructure(args.Get(mrna::bridge::OPT_RNASTRUCTURE_DATA), false);
#endif  // USE_RNASTRUCTURE

  const auto em = mrna::energy::EnergyModel::FromArgParse(args);

  auto cfg = mrna::fuzz::FuzzCfg::FromArgParse(args);
  const bool afl_mode = args.Has(mrna::OPT_AFL);

  if (afl_mode) {
#ifdef __AFL_FUZZ_TESTCASE_LEN
    __AFL_INIT();
    // This must be after __AFL_INIT and before __AFL_LOOP.
    auto buf = reinterpret_cast<const char*>(__AFL_FUZZ_TESTCASE_BUF);
    while (__AFL_LOOP(1000)) {
      int len = __AFL_FUZZ_TESTCASE_LEN;
      std::string data(buf, len);
      if (data.size() > 0) {
        cfg.seed = eng();
        mrna::Fuzzer fuzzer(mrna::StringToPrimary(data), cfg, em, rnastructure);
        const auto res = fuzzer.Run();
        if (!res.empty()) abort();
      }
    }
#endif
  } else {
    verify(args.PosSize() == 2, "require min and max length");
    const int min_len = args.Pos<int>(0);
    const int max_len = args.Pos<int>(1);
    const auto interval = args.Get<int>(OPT_PRINT_INTERVAL);

    verify(min_len > 0, "invalid min length");
    verify(max_len >= min_len, "invalid max len");
    std::uniform_int_distribution<int> len_dist(min_len, max_len);

    printf("Fuzzing [%d, %d] len RNAs - %s\n", min_len, max_len, cfg.Desc().c_str());

    // Normal mode.
    auto start_time = std::chrono::steady_clock::now();
    for (int64_t i = 0;; ++i) {
      if (interval > 0 &&
          std::chrono::duration_cast<std::chrono::seconds>(
              std::chrono::steady_clock::now() - start_time)
                  .count() > interval) {
        printf("Fuzzed %" PRId64 " RNA\n", i);
        start_time = std::chrono::steady_clock::now();
      }
      int len = len_dist(eng);
      auto r = mrna::Primary::Random(len);

      // TODO: respect OPT_RANDOM
      mrna::fuzz::Fuzzer fuzzer(std::move(r), cfg, em);
      const auto res = fuzzer.Run();
      if (!res.empty()) {
        for (const auto& s : res) printf("%s\n", s.c_str());
        printf("\n");
      }
    }
  }
}
