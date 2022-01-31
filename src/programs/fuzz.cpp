// Copyright 2016 Eliot Courtney.
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
    mrna::Opt().LongName("print-interval").Default("5").Help("status update every n seconds");

// Energy model options:
// TODO: Get this to work for RNAstructure.
inline const auto OPT_RANDOM = mrna::Opt().LongName("random").Help("use random energy models");

class FuzzHarness {
 public:
  explicit FuzzHarness(mrna::ArgParse args)
      : args_(std::move(args)), em_(mrna::energy::EnergyModel::FromArgParse(args_)),
        cfg_(mrna::fuzz::FuzzCfg::FromArgParse(args_)), e_(uint_fast32_t(time(nullptr))) {
#ifdef USE_RNASTRUCTURE
    rnastructure_ = std::make_unique<mrna::bridge::RNAstructure>(
        args.Get(mrna::bridge::OPT_RNASTRUCTURE_DATA), false);
#endif  // USE_RNASTRUCTURE
  }

  void Run() {
    printf("Fuzzing with config: %s\n", cfg_.Desc().c_str());

    if (args_.Has(mrna::OPT_AFL)) {
      DoAflFuzz();
    } else {
      DoRegularFuzz();
    }
  }

 private:
  mrna::ArgParse args_;
  mrna::energy::EnergyModel em_;
  mrna::fuzz::FuzzCfg cfg_;
  std::mt19937 e_;
  std::unique_ptr<mrna::bridge::RnaPackage> rnastructure_;

  void DoAflFuzz() {
#ifdef __AFL_FUZZ_TESTCASE_LEN
    __AFL_INIT();
    // This must be after __AFL_INIT and before __AFL_LOOP.
    auto buf = reinterpret_cast<const char*>(__AFL_FUZZ_TESTCASE_BUF);
    while (__AFL_LOOP(1000)) {
      int len = __AFL_FUZZ_TESTCASE_LEN;
      std::string data(buf, len);
      if (data.size() > 0) {
        cfg.seed = e_();
        mrna::Fuzzer fuzzer(mrna::StringToPrimary(data), em_, cfg_);
        const auto res = fuzzer.Run();
        if (!res.empty()) abort();
      }
    }
#endif
  }

  void DoRegularFuzz() {
    verify(args_.PosSize() == 2, "require min and max length");
    const int min_len = args_.Pos<int>(0);
    const int max_len = args_.Pos<int>(1);
    const auto interval = args_.Get<int>(OPT_PRINT_INTERVAL);

    verify(min_len > 0, "invalid min length");
    verify(max_len >= min_len, "invalid max len");
    std::uniform_int_distribution<int> len_dist(min_len, max_len);

    printf("Regular fuzzing [%d, %d] len RNAs\n", min_len, max_len);

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
      int len = len_dist(e_);
      auto r = mrna::Primary::Random(len);

      // TODO: respect OPT_RANDOM
      mrna::fuzz::Fuzzer fuzzer(std::move(r), em_, cfg_);
      const auto res = fuzzer.Run();
      if (!res.empty()) {
        for (const auto& s : res) printf("%s\n", s.c_str());
        printf("\n");
      }
    }
  }
};

int main(int argc, char* argv[]) {
  mrna::ArgParse args;
  mrna::fuzz::RegisterOpts(&args);
  args.RegisterOpt(mrna::OPT_AFL);
  args.RegisterOpt(mrna::bridge::OPT_RNASTRUCTURE_DATA);
  args.RegisterOpt(OPT_PRINT_INTERVAL);
  args.RegisterOpt(OPT_RANDOM);

  args.ParseOrExit(argc, argv);

  auto runner = FuzzHarness(std::move(args));
  runner.Run();
}
