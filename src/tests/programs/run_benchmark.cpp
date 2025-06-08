// Copyright 2025 Eliot Courtney.
#include <benchmark/benchmark.h>

#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/variadic/to_seq.hpp>
#include <random>
#include <tuple>
#include <utility>

#include "api/ctx/ctx.h"
#include "api/ctx/ctx_cfg.h"
#include "api/subopt/subopt_cfg.h"
#include "model/energy.h"
#include "model/primary.h"
#include "tests/init.h"
#include "util/argparse.h"
#include "util/error.h"

namespace mrna {

template <class... Args>
void Mfe(benchmark::State& state, Args&&... arglist) {
  auto args = std::make_tuple(std::move(arglist)...);
  const Ctx ctx(*std::get<0>(args), CtxCfg{.mfe_alg = std::get<1>(args)});
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
    auto result = ctx.SuboptIntoVector(r, cfg);
    benchmark::DoNotOptimize(result);
    benchmark::ClobberMemory();
  }
  state.SetComplexityN(state.range(0));
}

template <class... Args>
void Pfn(benchmark::State& state, Args&&... arglist) {
  auto args = std::make_tuple(std::move(arglist)...);
  const Ctx ctx(*std::get<0>(args), CtxCfg{.pfn_alg = std::get<1>(args)});
  std::mt19937 eng(0);

  for (auto _ : state) {
    auto r = Primary::Random(static_cast<int>(state.range(0)), eng);
    auto result = ctx.Pfn(r);
    benchmark::DoNotOptimize(result);
    benchmark::ClobberMemory();
  }
  state.SetComplexityN(state.range(0));
}

#define BENCHMARK_NAME(m, kind, ...) BOOST_PP_SEQ_CAT((m)(_)(kind))

#define FORCE_EVAL(...) __VA_ARGS__

#define DEFINE_MFE_BENCH1(r, m, kind)                        \
  BENCHMARK_CAPTURE(Mfe, m kind, &(m), CtxCfg::MfeAlg::kind) \
      ->RangeMultiplier(2)                                   \
      ->Range(16, 512)                                       \
      ->Complexity()                                         \
      ->Unit(benchmark::kMillisecond);

#define DEFINE_MFE_BENCH(m, ...) \
  BOOST_PP_SEQ_FOR_EACH(DEFINE_MFE_BENCH1, m, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))

#define DEFINE_SUBOPT_BENCH1(r, m, kind)                                                         \
  BENCHMARK_CAPTURE(                                                                             \
      Subopt, m kind 100strucs, &(m), CtxCfg::SuboptAlg::kind, subopt::SuboptCfg{.strucs = 100}) \
      ->RangeMultiplier(2)                                                                       \
      ->Range(16, 512)                                                                           \
      ->Complexity()                                                                             \
      ->Unit(benchmark::kMillisecond);                                                           \
  BENCHMARK_CAPTURE(                                                                             \
      Subopt, m kind delta, &(m), CtxCfg::SuboptAlg::kind, subopt::SuboptCfg{.delta = E(0.2)})   \
      ->RangeMultiplier(2)                                                                       \
      ->Range(16, 512)                                                                           \
      ->Complexity()                                                                             \
      ->Unit(benchmark::kMillisecond);

#define DEFINE_SUBOPT_BENCH(m, ...) \
  BOOST_PP_SEQ_FOR_EACH(DEFINE_SUBOPT_BENCH1, m, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))

#define DEFINE_PARTITION_BENCH1(r, m, kind)                  \
  BENCHMARK_CAPTURE(Pfn, m kind, &(m), CtxCfg::PfnAlg::kind) \
      ->RangeMultiplier(2)                                   \
      ->Range(16, 512)                                       \
      ->Complexity()                                         \
      ->Unit(benchmark::kMillisecond);

#define DEFINE_PARTITION_BENCH(m, ...) \
  BOOST_PP_SEQ_FOR_EACH(DEFINE_PARTITION_BENCH1, m, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))

DEFINE_MFE_BENCH(base_t04, DEBUG, SPARSE_OPT);
DEFINE_SUBOPT_BENCH(base_t04, DEBUG, ITERATIVE, ITERATIVE_LOWMEM, PERSISTENT, PERSISTENT_LOWMEM);
DEFINE_PARTITION_BENCH(base_t04, DEBUG, OPT);

DEFINE_MFE_BENCH(baseopt_t04, DEBUG, SPARSE_OPT);
DEFINE_SUBOPT_BENCH(baseopt_t04, DEBUG, ITERATIVE, ITERATIVE_LOWMEM, PERSISTENT, PERSISTENT_LOWMEM);
DEFINE_PARTITION_BENCH(baseopt_t04, DEBUG, OPT);

DEFINE_MFE_BENCH(stack_t04, DEBUG);
DEFINE_SUBOPT_BENCH(stack_t04, ITERATIVE, ITERATIVE_LOWMEM, PERSISTENT, PERSISTENT_LOWMEM);
// TODO(2): Add when partition is implemented for stack model.

// NEWMODEL/NEWBACKEND: Add new benches here.

}  // namespace mrna

int main(int argc, char** argv) {
  mrna::InitProgram();
  ::benchmark::Initialize(&argc, argv);
  mrna::ArgParse args;
  mrna::RegisterOptsBackendCfg(&args);
  args.ParseOrExit(argc, argv);
  mrna::InitTest(args.Get(mrna::OPT_MEMERNA_DATA));

  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();
}
