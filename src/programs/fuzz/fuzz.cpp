// Copyright 2016 Eliot Courtney.
#include <fmt/core.h>

#include <chrono>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "api/energy/pseudofree_cfg.h"
#include "fuzz/fuzz_cfg.h"
#include "fuzz/fuzz_harness.h"
#include "fuzz/fuzz_invocation.h"
#include "model/energy.h"
#include "model/primary.h"
#include "programs/print.h"
#include "util/argparse.h"
#include "util/error.h"
#include "util/log.h"

inline const auto OPT_PRINT_INTERVAL = mrna::Opt(mrna::Opt::ARG)
                                           .LongName("print-interval")
                                           .Default("5")
                                           .Help("status update every n seconds");
inline const auto OPT_ENUMERATE =
    mrna::Opt(mrna::Opt::FLAG).LongName("enumerate").Help("enumerate all sequences");

class FuzzRunner {
 public:
  explicit FuzzRunner(mrna::fuzz::FuzzCfg cfg, int interval)
      : cfg_(std::move(cfg)), harness_(cfg_), start_time_(Clock::now()),
        last_status_time_(start_time_), interval_(interval) {}

  bool RunSingle(const mrna::Primary& r) {
    mrna::loginfo("Running single fuzz on {}", r.ToSeq());
    for (int64_t i = 0;; ++i) {
      if (!KeepRunning(i, "times")) return true;
      if (RunInvocation(r)) return false;
    }
  }

  bool RunEnumerate(int min_len, int max_len) {
    mrna::loginfo("Exhaustive fuzzing [{}-{}] len RNAs", min_len, max_len);
    int64_t i = 0;
    mrna::Primary r(min_len);
    while (static_cast<int>(r.size()) <= max_len) {
      if (!KeepRunning(i)) return true;
      if (++i % 1000 == 0) mrna::loginfo("Fuzzed {} RNA, current size: {}", i, r.size());
      if (RunInvocation(r)) return false;
      r.Increment();
    }
    mrna::loginfo("Finished exhaustive fuzzing [{}-{}] len RNAs", min_len, max_len);
    return true;
  }

  bool RunRandom(int min_len, int max_len) {
    mrna::loginfo("Random fuzzing [{}-{}] len RNAs", min_len, max_len);
    std::uniform_int_distribution<int> len_dist(min_len, max_len);
    for (int64_t i = 0;; ++i) {
      if (!KeepRunning(i, "RNAs")) return true;
      const int len = len_dist(harness_.e());
      if (RunInvocation(mrna::Primary::Random(len, harness_.e()))) return false;
    }
  }

 private:
  using Clock = std::chrono::steady_clock;

  mrna::fuzz::FuzzCfg cfg_;
  mrna::fuzz::FuzzHarness harness_;
  Clock::time_point start_time_;
  Clock::time_point last_status_time_;
  int interval_;

  bool RunInvocation(const mrna::Primary& r) {
    mrna::erg::PseudofreeCfg pf(
        MaybeGetPairedPseudofree(r.size()), MaybeGetUnpairedPseudofree(r.size()));
    const auto res = harness_.Run(r, pf);
    MaybePrintResult(res, pf);
    return !res.empty();
  }

  std::vector<mrna::Energy> MaybeGetPairedPseudofree(std::size_t length) {
    if (!cfg_.pf_paired.empty()) return cfg_.pf_paired;
    return MaybeGetPseudofree(length);
  }

  std::vector<mrna::Energy> MaybeGetUnpairedPseudofree(std::size_t length) {
    if (!cfg_.pf_unpaired.empty()) return cfg_.pf_unpaired;
    return MaybeGetPseudofree(length);
  }

  std::vector<mrna::Energy> MaybeGetPseudofree(std::size_t length) {
    if (!cfg_.random_pseudofree) return {};
    return mrna::RandomEnergies(length, cfg_.random_pseudofree_cfg.min_energy,
        cfg_.random_pseudofree_cfg.max_energy, harness_.e());
  }

  [[nodiscard]] bool KeepRunning(int64_t i, const char* unit = nullptr) {
    const auto now = Clock::now();
    if (cfg_.fuzz_time_secs.has_value() &&
        std::chrono::duration_cast<std::chrono::seconds>(now - start_time_).count() >=
            *cfg_.fuzz_time_secs)
      return false;
    if (unit != nullptr && interval_ > 0 &&
        std::chrono::duration_cast<std::chrono::seconds>(now - last_status_time_).count() >
            interval_) {
      mrna::loginfo("Fuzzed {} {}", i, unit);
      last_status_time_ = now;
    }
    return true;
  }

  void MaybePrintResult(const mrna::fuzz::Error& res, const mrna::erg::PseudofreeCfg& pf) {
    if (res.empty()) return;
    fmt::print("Energy model: {}\n", cfg_.energy_model);
    fmt::print("Energy cfg: {}\n", cfg_.energy_cfg);
    fmt::print("Backends:{}\n", FormatBackends());
    if (auto random_cfg = harness_.last_random_model_cfg())
      fmt::print("Random model cfg: {}\n", *random_cfg);
    if (cfg_.random_pseudofree)
      fmt::print("Random pseudofree cfg: {}\n", cfg_.random_pseudofree_cfg);
    if (!pf.paired.empty())
      fmt::print("Pseudofree paired energies: {}\n", mrna::FormatPseudofreeEnergies(pf.paired));
    if (!pf.unpaired.empty())
      fmt::print("Pseudofree unpaired energies: {}\n", mrna::FormatPseudofreeEnergies(pf.unpaired));
    for (const auto& s : res) fmt::print("{}\n", s);
    fmt::print("\n");
  }

  [[nodiscard]] std::string FormatBackends() const {
    std::string backends;
    for (const auto& backend : cfg_.backends) backends += fmt::format(" {}", backend);
    return backends;
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
  const auto enumerate = args.GetOr(OPT_ENUMERATE);
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
  auto runner = FuzzRunner(fuzz_cfg, interval);
  bool ok = false;
  if (!seq.empty()) {
    ok = runner.RunSingle(mrna::Primary::FromSeq(seq));
  } else if (enumerate) {
    ok = runner.RunEnumerate(min_len, max_len);
  } else {
    ok = runner.RunRandom(min_len, max_len);
  }
  return ok ? 0 : 1;
}
