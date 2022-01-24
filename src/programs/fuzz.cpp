// Copyright 2016 E.
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <random>
#include <set>

#include "bridge/bridge.h"
#include "bridge/memerna.h"
#include "compute/energy/structure.h"
#include "compute/mfe/brute.h"
#include "compute/partition/brute.h"
#include "compute/subopt/brute.h"
#include "fuzz/fuzzer.h"
#include "util/string.h"

#ifdef USE_RNASTRUCTURE
#include "bridge/rnastructure.h"
#endif  // USE_RNASTRUCTURE

using mrna::Opt;

#ifdef __AFL_FUZZ_TESTCASE_LEN
#include <unistd.h>  // For __AFL_FUZZ_TESTCASE_LEN
__AFL_FUZZ_INIT();
#endif

int main(int argc, char* argv[]) {
  std::mt19937 eng(uint_fast32_t(time(nullptr)));
  mrna::ArgParse args({
      {"print-interval", Opt("status update every n seconds").Arg("-1")},
      {"afl", Opt("reads one rna from stdin and fuzzes - useful for use with afl")},
  });
  args.AddOptions(mrna::fuzz::FUZZ_OPTS);
  args.AddOptions(mrna::energy::ENERGY_OPTS);
  args.ParseOrExit(argc, argv);

  // TODO: Actually set this.
#ifdef USE_RNASTRUCTURE
  mrna::bridge::RNAstructure rnastructure(args.GetOption("rnastructure-data"), false);
#endif  // USE_RNASTRUCTURE

  const auto em = mrna::energy::EnergyModel::FromArgParse(args);

  auto cfg = mrna::fuzz::FuzzCfg::FromArgParse(args);
  const bool afl_mode = args.HasFlag("afl");
  verify(!cfg.mfe_rnastructure || !args.HasFlag("seed"),
      "seed option incompatible with rnastructure testing");

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
    auto pos = args.GetPositional();
    verify(pos.size() == 2, "require min and max length");
    const int min_len = atoi(pos[0].c_str());
    const int max_len = atoi(pos[1].c_str());
    const auto interval = atoi(args.GetOption("print-interval").c_str());

    verify(min_len > 0, "invalid min length");
    verify(max_len >= min_len, "invalid max len");
    std::uniform_int_distribution<int> len_dist(min_len, max_len);

    printf("Fuzzing [%d, %d] len RNAs - %s\n", min_len, max_len, cfg.Describe().c_str());

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

      cfg.seed = eng();
      mrna::fuzz::Fuzzer fuzzer(std::move(r), cfg, em);
      const auto res = fuzzer.Run();
      if (!res.empty()) {
        for (const auto& s : res) printf("%s\n", s.c_str());
        printf("\n");
      }
    }
  }
}
