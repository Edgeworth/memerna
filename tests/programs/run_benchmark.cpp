// Copyright 2022 Eliot Courtney.
#include <benchmark/benchmark.h>

#include <ios>

#include "common_test.h"
#include "compute/subopt/subopt_cfg.h"
#include "ctx/ctx.h"
#include "ctx/ctx_cfg.h"
#include "util/argparse.h"
#include "util/error.h"

namespace mrna {

template <class... Args>
void Mfe(benchmark::State& state, Args&&... arglist) {
  auto args = std::make_tuple(std::move(arglist)...);
  ctx::Ctx ctx(*std::get<0>(args), ctx::CtxCfg{.dp_alg = std::get<1>(args)});

  for (auto _ : state) {
    auto r = Primary::Random(static_cast<int>(state.range(0)));
    benchmark::DoNotOptimize(ctx.Fold(r));
    benchmark::ClobberMemory();
  }
  state.SetComplexityN(state.range(0));
}

template <class... Args>
void Subopt(benchmark::State& state, Args&&... arglist) {
  auto args = std::make_tuple(std::move(arglist)...);
  ctx::Ctx ctx(*std::get<0>(args), ctx::CtxCfg{.subopt_alg = std::get<1>(args)});
  auto cfg = std::get<2>(args);

  for (auto _ : state) {
    auto r = Primary::Random(static_cast<int>(state.range(0)));
    benchmark::DoNotOptimize(ctx.SuboptimalIntoVector(r, cfg));
    benchmark::ClobberMemory();
  }
  state.SetComplexityN(state.range(0));
}

template <class... Args>
void Partition(benchmark::State& state, Args&&... arglist) {
  auto args = std::make_tuple(std::move(arglist)...);
  ctx::Ctx ctx(*std::get<0>(args), ctx::CtxCfg{.part_alg = std::get<1>(args)});

  for (auto _ : state) {
    auto r = Primary::Random(static_cast<int>(state.range(0)));
    benchmark::DoNotOptimize(ctx.Partition(r));
    benchmark::ClobberMemory();
  }
  state.SetComplexityN(state.range(0));
}

#define DEFINE_BENCHES(em)                                                                  \
  BENCHMARK_CAPTURE(Mfe, em##_slowest, &(em), ctx::CtxCfg::DpAlg::SLOWEST)                  \
      ->RangeMultiplier(2)                                                                  \
      ->Range(16, 512)                                                                      \
      ->Complexity()                                                                        \
      ->Unit(benchmark::kMillisecond);                                                      \
                                                                                            \
  BENCHMARK_CAPTURE(Mfe, em##_fastest, &(em), ctx::CtxCfg::DpAlg::FASTEST)                  \
      ->RangeMultiplier(2)                                                                  \
      ->Range(16, 512)                                                                      \
      ->Complexity()                                                                        \
      ->Unit(benchmark::kMillisecond);                                                      \
                                                                                            \
  BENCHMARK_CAPTURE(Subopt, em##_slowest_100strucs, &(em), ctx::CtxCfg::SuboptAlg::SLOWEST, \
      subopt::SuboptCfg{.strucs = 100})                                                     \
      ->RangeMultiplier(2)                                                                  \
      ->Range(16, 512)                                                                      \
      ->Complexity()                                                                        \
      ->Unit(benchmark::kMillisecond);                                                      \
  BENCHMARK_CAPTURE(Subopt, em##_slowest_delta, &(em), ctx::CtxCfg::SuboptAlg::SLOWEST,     \
      subopt::SuboptCfg{.delta = E(0.2)})                                                   \
      ->RangeMultiplier(2)                                                                  \
      ->Range(16, 512)                                                                      \
      ->Complexity()                                                                        \
      ->Unit(benchmark::kMillisecond);                                                      \
                                                                                            \
  BENCHMARK_CAPTURE(Subopt, em##_fastest_100strucs, &(em), ctx::CtxCfg::SuboptAlg::FASTEST, \
      subopt::SuboptCfg{.strucs = 100})                                                     \
      ->RangeMultiplier(2)                                                                  \
      ->Range(16, 512)                                                                      \
      ->Complexity()                                                                        \
      ->Unit(benchmark::kMillisecond);                                                      \
                                                                                            \
  BENCHMARK_CAPTURE(Subopt, em##_fastest_delta, &(em), ctx::CtxCfg::SuboptAlg::FASTEST,     \
      subopt::SuboptCfg{.delta = E(0.2)})                                                   \
      ->RangeMultiplier(2)                                                                  \
      ->Range(16, 512)                                                                      \
      ->Complexity()                                                                        \
      ->Unit(benchmark::kMillisecond);                                                      \
                                                                                            \
  BENCHMARK_CAPTURE(Partition, em##_slowest, &(em), ctx::CtxCfg::PartAlg::SLOWEST)          \
      ->RangeMultiplier(2)                                                                  \
      ->Range(16, 512)                                                                      \
      ->Complexity()                                                                        \
      ->Unit(benchmark::kMillisecond);                                                      \
                                                                                            \
  BENCHMARK_CAPTURE(Partition, em##_fastest, &(em), ctx::CtxCfg::PartAlg::FASTEST)          \
      ->RangeMultiplier(2)                                                                  \
      ->Range(16, 512)                                                                      \
      ->Complexity()                                                                        \
      ->Unit(benchmark::kMillisecond)

#if ENERGY_PRECISION == 1

DEFINE_BENCHES(t04p1);

#elif ENERGY_PRECISION == 2

DEFINE_BENCHES(t04p2);

// NEWMODEL: Add new benches here.

#endif

}  // namespace mrna

int main(int argc, char** argv) {
  std::ios_base::sync_with_stdio(false);
  ::benchmark::Initialize(&argc, argv);
  mrna::ArgParse args;
  mrna::erg::RegisterOpts(&args);
  args.ParseOrExit(argc, argv);
  mrna::InitTest(args.Get(mrna::erg::OPT_MEMERNA_DATA));

  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();
}
