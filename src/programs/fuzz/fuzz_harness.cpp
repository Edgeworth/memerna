// Copyright 2022 Eliot Courtney.
#include "programs/fuzz/fuzz_harness.h"

#include <cinttypes>
#include <ctime>
#include <string>
#include <utility>

FuzzHarness::FuzzHarness(mrna::ArgParse args)
    : args_(std::move(args)), em_(mrna::energy::EnergyModel::FromArgParse(args_)),
      cfg_(mrna::fuzz::FuzzCfg::FromArgParse(args_)), e_(uint_fast32_t(time(nullptr))) {
#ifdef USE_RNASTRUCTURE
  rstr_ = std::make_shared<mrna::bridge::RNAstructure>(
      args_.Get(mrna::bridge::OPT_RNASTRUCTURE_DATA), false);
#endif  // USE_RNASTRUCTURE
  std::cout << "Fuzzing with config: " << cfg_.Desc() << '\n';
}

mrna::fuzz::FuzzInvocation FuzzHarness::CreateInvocation(const mrna::Primary& r) {
  mrna::fuzz::FuzzInvocation invoc(r, em_, cfg_);
  if (args_.GetOr(OPT_RANDOM_MODEL))
    invoc = mrna::fuzz::FuzzInvocation(r, mrna::energy::EnergyModel::Random(e_()), cfg_);
#ifdef USE_RNASTRUCTURE
  invoc.set_rnastructure(rstr_);
#endif  // USE_RNASTRUCTURE
  return invoc;
}
