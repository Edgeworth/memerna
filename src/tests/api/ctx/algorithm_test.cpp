// Copyright 2026 Eliot Courtney.
#include "api/ctx/algorithm.h"

#include <optional>
#include <string>
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
      /*kind=*/std::nullopt, /*alg=*/std::nullopt, kT04Cfg, erg::EnergyCfg{}, erg::PseudofreeCfg{});
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::BASE);
  EXPECT_EQ(res->alg, MfeAlg::SPARSE_OPT);
}

TEST(AlgorithmResolveTest, ResolveMfeAutoKeepsBaseWithPseudofree) {
  const auto res = ResolveMfe(
      /*kind=*/std::nullopt, /*alg=*/std::nullopt, kT04Cfg, erg::EnergyCfg{}, NonEmptyPf());
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::BASE);
  EXPECT_EQ(res->alg, MfeAlg::SPARSE_OPT);
}

TEST(AlgorithmResolveTest, ResolveMfeAutoSelectsStackOptForT22) {
  const auto res = ResolveMfe(
      /*kind=*/std::nullopt, /*alg=*/std::nullopt, kT22Cfg, erg::EnergyCfg{}, erg::PseudofreeCfg{});
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::STACK);
  EXPECT_EQ(res->alg, MfeAlg::OPT);
}

TEST(AlgorithmResolveTest, ResolveMfeCanonOptSupportsD2) {
  erg::EnergyCfg cfg;
  cfg.ctd = erg::EnergyCfg::Ctd::D2;
  const auto res = ResolveMfe(BackendKind::CANON, MfeAlg::OPT, kT04Cfg, cfg, erg::PseudofreeCfg{});
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::CANON);
  EXPECT_EQ(res->alg, MfeAlg::OPT);
}

TEST(AlgorithmResolveTest, ResolveMfeCanonOptRejectsT22) {
  const auto res =
      ResolveMfe(BackendKind::CANON, MfeAlg::OPT, kT22Cfg, erg::EnergyCfg{}, erg::PseudofreeCfg{});
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

TEST(AlgorithmResolveTest, ResolvePfnAutoSelectsBasePfnOptWhenSupported) {
  const auto res = ResolvePfn(
      /*kind=*/std::nullopt, /*alg=*/std::nullopt, kT04Cfg, erg::EnergyCfgPfn{},
      erg::PseudofreeCfg{});
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::BASE);
  EXPECT_EQ(res->alg, PfnAlg::OPT);
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
