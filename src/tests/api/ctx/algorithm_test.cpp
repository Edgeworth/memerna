// Copyright 2026 Eliot Courtney.
#include "api/ctx/algorithm.h"

#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "api/ctx/ctx.h"
#include "api/subopt/subopt_cfg.h"
#include "gtest/gtest.h"
#include "model/energy.h"
#include "model/primary.h"
#include "tests/init.h"

namespace mrna {

namespace {

erg::PseudofreeCfg NonEmptyPf() {
  return erg::PseudofreeCfg(std::vector<Energy>{E(1.0)}, std::vector<Energy>{E(1.0)});
}

bool HasBackendModel(const BackendModelPtr& m) {
  return std::visit([](const auto& model) { return model != nullptr; }, m);
}

}  // namespace

TEST(AlgorithmResolveTest, ResolveEfnAutoSelectsBaseWhenSupported) {
  const auto res =
      ResolveEfn(/*kind=*/std::nullopt, kT04Cfg, erg::EnergyCfg{}, erg::PseudofreeCfg{});
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::BASE);
}

TEST(AlgorithmResolveTest, ResolveEfnAutoKeepsBaseWithPseudofree) {
  const auto res = ResolveEfn(/*kind=*/std::nullopt, kT04Cfg, erg::EnergyCfg{}, NonEmptyPf());
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::BASE);
}

TEST(AlgorithmResolveTest, ResolveEfnAutoSelectsStackForT22) {
  const auto res =
      ResolveEfn(/*kind=*/std::nullopt, kT22Cfg, erg::EnergyCfg{}, erg::PseudofreeCfg{});
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::STACK);
}

TEST(AlgorithmResolveTest, ResolveEfnAutoFailsWhenNoBackendSupportsCfg) {
  erg::EnergyCfg cfg;
  cfg.ctd = erg::EnergyCfg::Ctd::D2;
  std::string log;
  const auto res =
      ResolveEfn(/*kind=*/std::nullopt, kT22Cfg, cfg, erg::PseudofreeCfg{}, /*log=*/&log);
  EXPECT_FALSE(res.has_value());
  EXPECT_FALSE(log.empty());
}

TEST(AlgorithmResolveTest, ResolveMfeAutoSelectsBaseSparseOptWhenSupported) {
  const auto res = ResolveMfe(
      /*kind=*/std::nullopt, /*alg=*/std::nullopt, kT04Cfg, erg::EnergyCfg{}, erg::PseudofreeCfg{},
      trace::TraceCfg{});
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::BASE);
  EXPECT_EQ(res->alg, MfeAlg::SPARSE_OPT);
}

TEST(AlgorithmResolveTest, ResolveMfeAutoKeepsBaseWithPseudofree) {
  const auto res = ResolveMfe(
      /*kind=*/std::nullopt, /*alg=*/std::nullopt, kT04Cfg, erg::EnergyCfg{}, NonEmptyPf(),
      trace::TraceCfg{});
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::BASE);
  EXPECT_EQ(res->alg, MfeAlg::SPARSE_OPT);
}

TEST(AlgorithmResolveTest, ResolveMfeAutoSelectsStackOptForT22) {
  const auto res = ResolveMfe(
      /*kind=*/std::nullopt, /*alg=*/std::nullopt, kT22Cfg, erg::EnergyCfg{}, erg::PseudofreeCfg{},
      trace::TraceCfg{});
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::STACK);
  EXPECT_EQ(res->alg, MfeAlg::OPT);
}

TEST(AlgorithmResolveTest, ResolveMfeCanonOptSupportsD2) {
  erg::EnergyCfg cfg;
  cfg.ctd = erg::EnergyCfg::Ctd::D2;
  const auto res = ResolveMfe(
      BackendKind::CANON, MfeAlg::OPT, kT04Cfg, cfg, erg::PseudofreeCfg{}, trace::TraceCfg{});
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::CANON);
  EXPECT_EQ(res->alg, MfeAlg::OPT);
}

TEST(AlgorithmResolveTest, ResolveMfeCanonOptRejectsT22) {
  const auto res = ResolveMfe(BackendKind::CANON, MfeAlg::OPT, kT22Cfg, erg::EnergyCfg{},
      erg::PseudofreeCfg{}, trace::TraceCfg{});
  EXPECT_FALSE(res.has_value());
}

TEST(AlgorithmResolveTest, ResolveSuboptAutoSelectsBaseIterativeWhenSupported) {
  const subopt::SuboptCfg subopt_cfg;
  const auto res = ResolveSubopt(/*kind=*/std::nullopt, /*alg=*/std::nullopt, kT04Cfg,
      erg::EnergyCfg{}, erg::PseudofreeCfg{}, subopt_cfg);
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::BASE);
  EXPECT_EQ(res->alg, SuboptAlg::ITERATIVE);
}

TEST(AlgorithmResolveTest, ResolveSuboptAutoKeepsBaseWithPseudofree) {
  const subopt::SuboptCfg subopt_cfg;
  const auto res = ResolveSubopt(/*kind=*/std::nullopt, /*alg=*/std::nullopt, kT04Cfg,
      erg::EnergyCfg{}, NonEmptyPf(), subopt_cfg);
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::BASE);
  EXPECT_EQ(res->alg, SuboptAlg::ITERATIVE);
}

TEST(AlgorithmResolveTest, ResolveSuboptCanonPersistentSupportsD2) {
  erg::EnergyCfg cfg;
  cfg.ctd = erg::EnergyCfg::Ctd::D2;
  const subopt::SuboptCfg subopt_cfg;
  const auto res = ResolveSubopt(
      BackendKind::CANON, SuboptAlg::PERSISTENT, kT04Cfg, cfg, erg::PseudofreeCfg{}, subopt_cfg);
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::CANON);
  EXPECT_EQ(res->alg, SuboptAlg::PERSISTENT);
}

TEST(CtxSuboptTest, CallbackCanStopEnumeration) {
  const Ctx ctx(base_t04, kT04Cfg);
  const auto r = Primary::FromSeq("CCUCCGGG");
  int64_t calls = 0;
  const auto count = ctx.Subopt(r, /*mfe_alg=*/std::nullopt, SuboptAlg::ITERATIVE, erg::EnergyCfg{},
      erg::PseudofreeCfg{}, subopt::SuboptCfg{}, [&](const subopt::SuboptResult&) {
        ++calls;
        return false;
      });
  EXPECT_EQ(count, 1);
  EXPECT_EQ(calls, 1);
}

TEST(CtxSuboptTest, BackendSelectionRejectsInvalidConfiguration) {
  const Ctx ctx(base_t04, kT04Cfg);
  const auto infinity = std::numeric_limits<float>::infinity();
  const auto nan = std::numeric_limits<float>::quiet_NaN();
  const std::vector<subopt::SuboptCfg> invalid_cfgs{
      subopt::SuboptCfg{.strucs = -1},
      subopt::SuboptCfg{.delta = E(-0.1)},
      subopt::SuboptCfg{.delta = Energy::FromRaw(CAP_E.v + 1)},
      subopt::SuboptCfg{.time_secs = infinity},
      subopt::SuboptCfg{.time_secs = -infinity},
      subopt::SuboptCfg{.time_secs = nan},
      subopt::SuboptCfg{.time_secs = -2.0F},
  };

  for (const auto& cfg : invalid_cfgs) {
    EXPECT_THROW((void)ctx.BackendForSubopt(SuboptAlg::DEBUG, /*mfe_alg=*/std::nullopt,
                     erg::EnergyCfg{}, erg::PseudofreeCfg{}, cfg),
        std::runtime_error);
  }
}

TEST(CtxFoldTest, RandomTraceAutoSelectsStack) {
  auto backend_cfg = kT04Cfg;
  backend_cfg.data_src = RandomModelCfg(0, E(-1.0), E(1.0));
  const Ctx ctx(/*backend=*/std::nullopt, backend_cfg);
  const trace::TraceCfg trace_cfg{.random = true, .seed = 0};

  const auto backend =
      ctx.BackendForFold(/*alg=*/std::nullopt, erg::EnergyCfg{}, erg::PseudofreeCfg{}, trace_cfg);
  EXPECT_EQ(GetBackendKind(backend.m), BackendKind::STACK);

  const auto r = Primary::FromSeq("CCUCCGGG");
  const auto result =
      ctx.Fold(r, /*alg=*/std::nullopt, erg::EnergyCfg{}, erg::PseudofreeCfg{}, trace_cfg);
  EXPECT_EQ(result.tb.s.size(), r.size());
}

TEST(CtxFoldTest, RandomTraceRejectsIncompatibleExplicitBackend) {
  const Ctx ctx(base_t04, kT04Cfg);
  EXPECT_THROW((void)ctx.BackendForFold(/*alg=*/std::nullopt, erg::EnergyCfg{},
                   erg::PseudofreeCfg{}, trace::TraceCfg{.random = true, .seed = std::nullopt}),
      std::runtime_error);
}

TEST(CtxFoldTest, BackendSelectionRejectsTraceSeedWithoutRandomTrace) {
  const Ctx ctx(base_t04, kT04Cfg);
  EXPECT_THROW((void)ctx.BackendForFold(/*alg=*/std::nullopt, erg::EnergyCfg{},
                   erg::PseudofreeCfg{}, trace::TraceCfg{.random = false, .seed = 0}),
      std::runtime_error);
}

TEST(CtxTest, EfnProducesSerializableCtdsForEnclosingPair) {
  const auto r = Primary::FromSeq("GAAAC");
  const auto s = Secondary::FromDb("(...)");

  for (const auto& model : {BackendModelPtr{base_t04}, BackendModelPtr{stack_t04}}) {
    constexpr erg::EnergyCfg cfg{};
    const auto result = TotalEnergy(
        model, r, s, /*given_ctd=*/nullptr, cfg, erg::PseudofreeCfg{}, /*build_structure=*/true);
    EXPECT_EQ("[...]", cfg.ToCtdString(s, result.ctd));
  }
}

TEST(CtxTest, RejectsEmptySequences) {
  const Primary r;
  const Secondary s;
  for (const auto& model : t04_ms) {
    const Ctx ctx(model, kT04Cfg);
    EXPECT_THROW((void)ctx.Efn(r, s, {}), std::runtime_error);
    EXPECT_THROW((void)ctx.Fold(r, std::nullopt, {}, {}, {}), std::runtime_error);
    EXPECT_THROW(
        (void)ctx.SuboptIntoVector(r, std::nullopt, std::nullopt, {}, {}, {}), std::runtime_error);
    EXPECT_THROW((void)ctx.Pfn(r, std::nullopt, {}, {}), std::runtime_error);
  }
}

TEST(CtxBackendTest, BackendDescriptorsOwnTheirModelVariant) {
  const auto fold = Ctx(base_t04, kT04Cfg)
                        .BackendForFold(/*alg=*/std::nullopt, erg::EnergyCfg{},
                            erg::PseudofreeCfg{}, trace::TraceCfg{});
  const auto subopt = Ctx(base_t04, kT04Cfg)
                          .BackendForSubopt(/*alg=*/std::nullopt, /*mfe_alg=*/std::nullopt,
                              erg::EnergyCfg{}, erg::PseudofreeCfg{}, subopt::SuboptCfg{});
  const auto pfn = Ctx(base_t04, kT04Cfg)
                       .BackendForPfn(
                           /*alg=*/std::nullopt, erg::EnergyCfgPfn{}, erg::PseudofreeCfg{});

  EXPECT_TRUE(HasBackendModel(fold.m));
  EXPECT_TRUE(HasBackendModel(subopt.m));
  EXPECT_TRUE(HasBackendModel(pfn.m));

  Ctx source(base_t04, kT04Cfg);
  const auto before_move = source.BackendForFold(
      /*alg=*/std::nullopt, erg::EnergyCfg{}, erg::PseudofreeCfg{}, trace::TraceCfg{});
  Ctx destination(std::move(source));
  EXPECT_TRUE(HasBackendModel(before_move.m));
  EXPECT_EQ(GetBackendKind(before_move.m), BackendKind::BASE);
}

TEST(AlgorithmResolveTest, ResolvePfnAutoSelectsBasePfnOptWhenSupported) {
  const auto res = ResolvePfn(
      /*kind=*/std::nullopt, /*alg=*/std::nullopt, kT04Cfg, erg::EnergyCfgPfn{},
      erg::PseudofreeCfg{});
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::BASE);
  EXPECT_EQ(res->alg, PfnAlg::OPT);
}

TEST(AlgorithmResolveTest, ResolvePfnExplicitlySelectsBrute) {
  const auto res = ResolvePfn(
      BackendKind::BASE, PfnAlg::BRUTE, kT04Cfg, erg::EnergyCfgPfn{}, erg::PseudofreeCfg{});
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::BASE);
  EXPECT_EQ(res->alg, PfnAlg::BRUTE);
}

TEST(AlgorithmResolveTest, ResolvePfnBruteRejectsT22) {
  const auto res = ResolvePfn(
      BackendKind::STACK, PfnAlg::BRUTE, kT22Cfg, erg::EnergyCfgPfn{}, erg::PseudofreeCfg{});
  EXPECT_FALSE(res.has_value());
}

TEST(AlgorithmResolveTest, ResolvePfnAutoKeepsBaseWithPseudofree) {
  const auto res = ResolvePfn(
      /*kind=*/std::nullopt, /*alg=*/std::nullopt, kT04Cfg, erg::EnergyCfgPfn{}, NonEmptyPf());
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::BASE);
  EXPECT_EQ(res->alg, PfnAlg::OPT);
}

TEST(AlgorithmResolveTest, ResolvePfnCanonOptSupportsD2) {
  erg::EnergyCfg cfg;
  cfg.ctd = erg::EnergyCfg::Ctd::D2;
  const auto res = ResolvePfn(
      BackendKind::CANON, PfnAlg::OPT, kT04Cfg, erg::EnergyCfgPfn{cfg}, erg::PseudofreeCfg{});
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::CANON);
  EXPECT_EQ(res->alg, PfnAlg::OPT);
}

TEST(AlgorithmResolveTest, ResolvePfnAutoFailsWhenNoBackendProvidesPfnAlg) {
  std::string log;
  const auto res = ResolvePfn(
      /*kind=*/std::nullopt, /*alg=*/std::nullopt, kT22Cfg, erg::EnergyCfgPfn{},
      erg::PseudofreeCfg{}, /*log=*/&log);
  EXPECT_FALSE(res.has_value());
  EXPECT_FALSE(log.empty());
}

}  // namespace mrna
