// Copyright 2022 E.
#include "programs/fuzz/fuzz_harness.h"

#include <fmt/core.h>

#include <cinttypes>
#include <ctime>
#include <string>
#include <utility>

#include "api/energy/energy_cfg.h"

FuzzHarness::FuzzHarness(mrna::ArgParse args)
    : args_(std::move(args)), em_(mrna::erg::FromArgParse(args_)),
      cfg_(mrna::fuzz::FuzzCfg::FromArgParse(args_)), e_(uint_fast32_t(time(nullptr))) {
#ifdef USE_RNASTRUCTURE
  rstr_ = std::make_shared<mrna::bridge::RNAstructure>(
      args_.Get(mrna::bridge::OPT_RNASTRUCTURE_DATA), false);
#endif  // USE_RNASTRUCTURE
  fmt::print("Fuzzing with config: {}\n", cfg_.Desc());
}

mrna::fuzz::FuzzInvocation FuzzHarness::CreateInvocation(const mrna::Primary& r) {
  mrna::fuzz::FuzzInvocation invoc(r, em_, cfg_, 0);
  if (args_.GetOr(OPT_RANDOM_MODEL)) {
    uint_fast32_t seed = e_();
    invoc = mrna::fuzz::FuzzInvocation(r,
        mrna::erg::Random(args_.Get<mrna::erg::ModelKind>(mrna::erg::OPT_ENERGY_MODEL), seed), cfg_,
        seed);
  }
#ifdef USE_RNASTRUCTURE
  invoc.set_rnastructure(rstr_);
#endif  // USE_RNASTRUCTURE
  return invoc;
}
