// Copyright 2016 Eliot Courtney.
#include <fmt/core.h>

#include <chrono>
#include <random>
#include <string>

#include "fuzz/fuzz_cfg.h"
#include "fuzz/fuzz_harness.h"
#include "fuzz/fuzz_invocation.h"
#include "model/energy.h"
#include "model/primary.h"
#include "util/argparse.h"
#include "util/error.h"

inline const auto OPT_PRINT_INTERVAL = mrna::Opt(mrna::Opt::ARG)
                                           .LongName("print-interval")
                                           .Default("5")
                                           .Help("status update every n seconds");
inline const auto OPT_ENUMERATE =
    mrna::Opt(mrna::Opt::FLAG).LongName("enumerate").Help("enumerate all sequences");

class FuzzRunner {
 public:
  explicit FuzzRunner(mrna::fuzz::FuzzCfg cfg) : cfg_(std::move(cfg)), harness_(cfg_) {}

  void RunSingle(const mrna::Primary& r) {
    fmt::print("Running single fuzz on {}\n", r.ToSeq());
    RunInvocation(r);
  }

  void RunEnumerate(int min_len, int max_len) {
    fmt::print("Exhaustive fuzzing [{}-{}] len RNAs\n", min_len, max_len);
    int64_t i = 0;
    mrna::Primary r(min_len);
    while (static_cast<int>(r.size()) <= max_len) {
      if (++i % 1000 == 0) fmt::print("Fuzzed {} RNA, current size: {}\n", i, r.size());
      RunInvocation(r);
      r.Increment();
    }
    fmt::print("Finished exhaustive fuzzing [{}-{}] len RNAs\n", min_len, max_len);
  }

  void RunRandom(int min_len, int max_len, int interval) {
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
      const int len = len_dist(harness_.e());
      RunInvocation(mrna::Primary::Random(len, harness_.e()));
    }
  }

 private:
  mrna::fuzz::FuzzCfg cfg_;
  mrna::fuzz::FuzzHarness harness_;

  void RunInvocation(const mrna::Primary& r) {
    auto pf_paired = MaybeGetPseudofree(r.size());
    auto pf_unpaired = MaybeGetPseudofree(r.size());
    auto invoc = harness_.CreateInvocation(r, pf_paired, pf_unpaired);
    MaybePrintResult(invoc.Run(), pf_paired, pf_unpaired);
  }

  std::vector<mrna::Energy> MaybeGetPseudofree(std::size_t length) {
    if (!cfg_.random_pseudofree) return {};
    return mrna::RandomEnergies(length, mrna::E(-100.0), mrna::E(100.0), harness_.e());
  }

  void MaybePrintResult(const mrna::fuzz::Error& res, const std::vector<mrna::Energy>& pf_paired,
      const std::vector<mrna::Energy>& pf_unpaired) {
    if (res.empty()) return;
    if (cfg_.random_models) fmt::print("Random model seed: {}\n", harness_.last_seed().value());
    if (!pf_paired.empty()) {
      fmt::print("Pseudofree paired energies: ");
      PrintPseudofreeEnergy(pf_paired);
    }
    if (!pf_unpaired.empty()) {
      fmt::print("Pseudofree unpaired energies: ");
      PrintPseudofreeEnergy(pf_unpaired);
    }
    for (const auto& s : res) fmt::print("{}\n", s);
    fmt::print("\n");
  }

  static void PrintPseudofreeEnergy(const std::vector<mrna::Energy>& pf) {
    bool first = true;
    for (const auto& e : pf) {
      if (!first) fmt::print(",");
      first = false;
      fmt::print("{}", e);
    }

    fmt::print("\n");
  }
};

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
  auto runner = FuzzRunner(fuzz_cfg);
  if (!seq.empty()) {
    runner.RunSingle(mrna::Primary::FromSeq(seq));
  } else if (enumerate) {
    runner.RunEnumerate(min_len, max_len);
  } else {
    runner.RunRandom(min_len, max_len, interval);
  }
}
