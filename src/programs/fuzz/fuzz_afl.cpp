// Copyright 2022 Eliot Courtney.
#include <utility>

#include "bridge/bridge.h"
#include "fuzz/config.h"
#include "programs/fuzz/fuzz_harness.h"
#include "util/argparse.h"

#ifdef __AFL_FUZZ_TESTCASE_LEN
#include <unistd.h>  // For __AFL_FUZZ_TESTCASE_LEN

__AFL_FUZZ_INIT();
#endif

int main(int argc, char* argv[]) {
  mrna::ArgParse args;
  mrna::fuzz::RegisterOpts(&args);
  args.RegisterOpt(mrna::bridge::OPT_RNASTRUCTURE_DATA);
  args.RegisterOpt(OPT_RANDOM);
  args.ParseOrExit(argc, argv);

  auto harness = FuzzHarness(std::move(args));

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
      try {
        auto invoc = harness.CreateInvocation(mrna::Primary::FromSeq(rs));
        const auto res = invoc.Run();
        if (!res.empty()) {
          for (const auto& s : res) printf("%s\n", s.c_str());
          printf("\n");
          abort();
        }
      } catch (const std::exception& e) {
        // Ignore. Probably a bad input.
      }
    }
  }
#endif
}
