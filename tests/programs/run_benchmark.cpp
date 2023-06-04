// Copyright 2022 Eliot Courtney.
#include <benchmark/benchmark.h>

#include <memory>
#include <random>
#include <tuple>
#include <utility>

#include "api/ctx/ctx.h"
#include "api/ctx/ctx_cfg.h"
#include "api/energy/energy.h"
#include "api/subopt/subopt_cfg.h"
#include "common_test.h"
#include "model/energy.h"
#include "model/primary.h"
#include "util/argparse.h"
#include "util/error.h"

namespace mrna {

template <class... Args>
void Mfe(benchmark::State& state, Args&&... arglist) {
  auto args = std::make_tuple(std::move(arglist)...);
  const Ctx ctx(*std::get<0>(args), CtxCfg{.dp_alg = std::get<1>(args)});
  std::mt19937 eng(0);

  for (auto _ : state) {
    auto r = Primary::Random(static_cast<int>(state.range(0)), eng);
    auto result = ctx.Fold(r, {});
    benchmark::DoNotOptimize(result);
    benchmark::ClobberMemory();
  }
  state.SetComplexityN(state.range(0));
}

template <class... Args>
void Subopt(benchmark::State& state, Args&&... arglist) {
  auto args = std::make_tuple(std::move(arglist)...);
  const Ctx ctx(*std::get<0>(args), CtxCfg{.subopt_alg = std::get<1>(args)});
  auto cfg = std::get<2>(args);
  std::mt19937 eng(0);

  for (auto _ : state) {
    auto r = Primary::Random(static_cast<int>(state.range(0)), eng);
    auto result = ctx.SuboptimalIntoVector(r, cfg);
    benchmark::DoNotOptimize(result);
    benchmark::ClobberMemory();
  }
  state.SetComplexityN(state.range(0));
}

template <class... Args>
void Partition(benchmark::State& state, Args&&... arglist) {
  auto args = std::make_tuple(std::move(arglist)...);
  const Ctx ctx(*std::get<0>(args), CtxCfg{.part_alg = std::get<1>(args)});
  std::mt19937 eng(0);

  for (auto _ : state) {
    auto r = Primary::Random(static_cast<int>(state.range(0)), eng);
    auto result = ctx.Partition(r);
    benchmark::DoNotOptimize(result);
    benchmark::ClobberMemory();
  }
  state.SetComplexityN(state.range(0));
}

#define DEFINE_BENCHES(em)                                                             \
  BENCHMARK_CAPTURE(Mfe, em##_slowest, &(em), CtxCfg::DpAlg::SLOWEST)                  \
      ->RangeMultiplier(2)                                                             \
      ->Range(16, 512)                                                                 \
      ->Complexity()                                                                   \
      ->Unit(benchmark::kMillisecond);                                                 \
                                                                                       \
  BENCHMARK_CAPTURE(Mfe, em##_fastest, &(em), CtxCfg::DpAlg::FASTEST)                  \
      ->RangeMultiplier(2)                                                             \
      ->Range(16, 512)                                                                 \
      ->Complexity()                                                                   \
      ->Unit(benchmark::kMillisecond);                                                 \
                                                                                       \
  BENCHMARK_CAPTURE(Subopt, em##_slowest_100strucs, &(em), CtxCfg::SuboptAlg::SLOWEST, \
      subopt::SuboptCfg{.strucs = 100})                                                \
      ->RangeMultiplier(2)                                                             \
      ->Range(16, 512)                                                                 \
      ->Complexity()                                                                   \
      ->Unit(benchmark::kMillisecond);                                                 \
  BENCHMARK_CAPTURE(Subopt, em##_slowest_delta, &(em), CtxCfg::SuboptAlg::SLOWEST,     \
      subopt::SuboptCfg{.delta = E(0.2)})                                              \
      ->RangeMultiplier(2)                                                             \
      ->Range(16, 512)                                                                 \
      ->Complexity()                                                                   \
      ->Unit(benchmark::kMillisecond);                                                 \
                                                                                       \
  BENCHMARK_CAPTURE(Subopt, em##_fastest_100strucs, &(em), CtxCfg::SuboptAlg::FASTEST, \
      subopt::SuboptCfg{.strucs = 100})                                                \
      ->RangeMultiplier(2)                                                             \
      ->Range(16, 512)                                                                 \
      ->Complexity()                                                                   \
      ->Unit(benchmark::kMillisecond);                                                 \
                                                                                       \
  BENCHMARK_CAPTURE(Subopt, em##_fastest_delta, &(em), CtxCfg::SuboptAlg::FASTEST,     \
      subopt::SuboptCfg{.delta = E(0.2)})                                              \
      ->RangeMultiplier(2)                                                             \
      ->Range(16, 512)                                                                 \
      ->Complexity()                                                                   \
      ->Unit(benchmark::kMillisecond);                                                 \
                                                                                       \
  BENCHMARK_CAPTURE(Partition, em##_slowest, &(em), CtxCfg::PartAlg::SLOWEST)          \
      ->RangeMultiplier(2)                                                             \
      ->Range(16, 512)                                                                 \
      ->Complexity()                                                                   \
      ->Unit(benchmark::kMillisecond);                                                 \
                                                                                       \
  BENCHMARK_CAPTURE(Partition, em##_fastest, &(em), CtxCfg::PartAlg::FASTEST)          \
      ->RangeMultiplier(2)                                                             \
      ->Range(16, 512)                                                                 \
      ->Complexity()                                                                   \
      ->Unit(benchmark::kMillisecond)

#if ENERGY_PRECISION == 1

DEFINE_BENCHES(t04p1);

#elif ENERGY_PRECISION == 2

DEFINE_BENCHES(t04p2);
DEFINE_BENCHES(t12p2);
// TODO(2): Uncomment when partition+subopt are implemented.
// DEFINE_BENCHES(t22p2);

// NEWMODEL: Add new benches here.

#endif

}  // namespace mrna

int main(int argc, char** argv) {
  mrna::InitProgram();
  ::benchmark::Initialize(&argc, argv);
  mrna::ArgParse args;
  mrna::erg::RegisterOptsEnergyModel(&args);
  args.ParseOrExit(argc, argv);
  mrna::InitTest(args.Get(mrna::erg::OPT_MEMERNA_DATA));

  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();
}
