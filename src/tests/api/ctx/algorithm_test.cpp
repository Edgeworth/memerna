// Copyright 2026 Eliot Courtney.
#include "api/ctx/algorithm.h"

#include <optional>
#include <string>
#include <vector>

#include "api/subopt/subopt_cfg.h"
#include "gtest/gtest.h"
#include "model/energy.h"
#include "tests/init.h"

namespace mrna {

namespace {

erg::PseudofreeCfg NonEmptyPf() {
  return erg::PseudofreeCfg(std::vector<Energy>{E(1.0)}, std::vector<Energy>{E(1.0)});
}

}  // namespace

TEST(AlgorithmResolveTest, ResolveEfnAutoPrefersBaseoptWhenSupported) {
  const auto res =
      ResolveEfn(/*kind=*/std::nullopt, kT04Cfg, erg::EnergyCfg{}, erg::PseudofreeCfg{});
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::BASEOPT);
}

TEST(AlgorithmResolveTest, ResolveEfnAutoFallsBackWhenBaseoptUnsupported) {
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

TEST(AlgorithmResolveTest, ResolveMfeAutoPrefersBaseoptSparseOptWhenSupported) {
  const auto res = ResolveMfe(
      /*kind=*/std::nullopt, /*alg=*/std::nullopt, kT04Cfg, erg::EnergyCfg{}, erg::PseudofreeCfg{});
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::BASEOPT);
  EXPECT_EQ(res->alg, MfeAlg::SPARSE_OPT);
}

TEST(AlgorithmResolveTest, ResolveMfeAutoFallsBackWhenBaseoptUnsupported) {
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

TEST(AlgorithmResolveTest, ResolveSuboptAutoPrefersBaseoptIterativeWhenSupported) {
  const subopt::SuboptCfg subopt_cfg;
  const auto res = ResolveSubopt(/*kind=*/std::nullopt, /*alg=*/std::nullopt, kT04Cfg,
      erg::EnergyCfg{}, erg::PseudofreeCfg{}, subopt_cfg);
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::BASEOPT);
  EXPECT_EQ(res->alg, SuboptAlg::ITERATIVE);
}

TEST(AlgorithmResolveTest, ResolveSuboptAutoFallsBackWhenBaseoptUnsupported) {
  const subopt::SuboptCfg subopt_cfg;
  const auto res = ResolveSubopt(/*kind=*/std::nullopt, /*alg=*/std::nullopt, kT04Cfg,
      erg::EnergyCfg{}, NonEmptyPf(), subopt_cfg);
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::BASE);
  EXPECT_EQ(res->alg, SuboptAlg::ITERATIVE);
}

TEST(AlgorithmResolveTest, ResolvePfnAutoPrefersBaseoptOptWhenSupported) {
  const auto res = ResolvePfn(
      /*kind=*/std::nullopt, /*alg=*/std::nullopt, kT04Cfg, erg::EnergyCfg{}, erg::PseudofreeCfg{});
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::BASEOPT);
  EXPECT_EQ(res->alg, PfnAlg::OPT);
}

TEST(AlgorithmResolveTest, ResolvePfnAutoFallsBackWhenBaseoptUnsupported) {
  const auto res = ResolvePfn(
      /*kind=*/std::nullopt, /*alg=*/std::nullopt, kT04Cfg, erg::EnergyCfg{}, NonEmptyPf());
  ASSERT_TRUE(res.has_value());
  EXPECT_EQ(res->backend, BackendKind::BASE);
  EXPECT_EQ(res->alg, PfnAlg::OPT);
}

TEST(AlgorithmResolveTest, ResolvePfnAutoFailsWhenNoBackendProvidesPfnAlg) {
  std::string log;
  const auto res = ResolvePfn(
      /*kind=*/std::nullopt, /*alg=*/std::nullopt, kT22Cfg, erg::EnergyCfg{}, erg::PseudofreeCfg{},
      /*log=*/&log);
  EXPECT_FALSE(res.has_value());
  EXPECT_FALSE(log.empty());
}

}  // namespace mrna
