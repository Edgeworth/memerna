// Copyright 2022 Eliot Courtney.
#ifndef PROGRAMS_FUZZ_FUZZ_HARNESS_H_
#define PROGRAMS_FUZZ_FUZZ_HARNESS_H_
#include <memory>
#include <random>

#include "api/energy/model.h"
#include "fuzz/fuzz_cfg.h"
#include "fuzz/fuzz_invocation.h"
#include "model/primary.h"
#include "util/argparse.h"

#ifdef USE_RNASTRUCTURE
#include "api/bridge/rnastructure.h"
#endif  // USE_RNASTRUCTURE

// Energy model options:
inline const auto OPT_RANDOM_MODELS =
    mrna::Opt(mrna::Opt::FLAG)
        .LongName("random-models")
        .Help("use random energy models - different each time. c.f. seed");

class FuzzHarness {
 public:
  explicit FuzzHarness(mrna::ArgParse args);

  mrna::fuzz::FuzzInvocation CreateInvocation(const mrna::Primary& r);
  std::mt19937& e() { return e_; }

 private:
  mrna::ArgParse args_;
  mrna::erg::EnergyModelPtr em_;
  mrna::fuzz::FuzzCfg cfg_;
  std::mt19937 e_;

#ifdef USE_RNASTRUCTURE
  std::shared_ptr<mrna::bridge::RNAstructure> rstr_;
#endif  // USE_RNASTRUCTURE
};

#endif  // PROGRAMS_FUZZ_FUZZ_HARNESS_H_
