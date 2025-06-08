// Copyright 2022 Eliot Courtney.
#include <string>
#include <utility>

#include "fuzz/fuzz_cfg.h"
#include "fuzz/fuzz_harness.h"
#include "util/argparse.h"
#include "util/error.h"

#ifdef __AFL_FUZZ_TESTCASE_LEN
#include <unistd.h>  // For __AFL_FUZZ_TESTCASE_LEN

__AFL_FUZZ_INIT();
#endif

inline const auto OPT_MAX_LEN = mrna::Opt(mrna::Opt::ARG)
                                    .LongName("max-len")
                                    .Help("limit max length of sequences fuzzed")
                                    .Default(-1);

int main(int argc, char* argv[]) {
  mrna::InitProgram();
  mrna::ArgParse args;
  mrna::fuzz::RegisterOpts(&args);
  args.RegisterOpt(OPT_MAX_LEN);
  args.ParseOrExit(argc, argv);

  [[maybe_unused]] const auto max_len = args.Get<int>(OPT_MAX_LEN);

  auto fuzz_cfg = mrna::fuzz::FuzzCfg::FromArgParse(args);
  auto harness = mrna::fuzz::FuzzHarness(std::move(fuzz_cfg));

#ifdef __AFL_FUZZ_TESTCASE_LEN
  __AFL_INIT();
  // This must be after __AFL_INIT and before __AFL_LOOP.
  auto buf = reinterpret_cast<const char*>(__AFL_FUZZ_TESTCASE_BUF);
  while (__AFL_LOOP(1000)) {
    int len = __AFL_FUZZ_TESTCASE_LEN;
    std::string data(buf, len);
    std::stringstream ss(data);
    std::string rs;
    while (ss >> rs) {
      mrna::Primary seq;
      try {
        if (max_len > 0 && static_cast<int>(rs.size()) > max_len) rs.resize(max_len);
        seq = mrna::Primary::FromSeq(rs);
      } catch (const std::exception& e) {
        // Ignore. Probably a bad input.
        continue;
      }
      auto invoc = harness.CreateInvocation(seq, {}, {});
      const auto res = invoc.Run();
      if (!res.empty()) {
        for (const auto& s : res) fmt::print("{}\n", s);
        fmt::print("\n");
        abort();
      }
    }
  }
#endif
}
