// Copyright 2024 Eliot Courtney.
#ifndef TESTS_UTIL_H_
#define TESTS_UTIL_H_

#include <string>
#include <tuple>

#include "api/ctx/backend.h"
#include "api/ctx/ctx.h"
#include "api/ctx/ctx_cfg.h"
#include "gtest/gtest.h"
#include "model/primary.h"
#include "model/secondary.h"

namespace mrna {

#define EXPECT_REL_EQ(a, b)                                                    \
  do {                                                                         \
    auto acopy = (a);                                                          \
    auto bcopy = (b);                                                          \
    EXPECT_TRUE(rel_eq(acopy, bcopy))                                          \
        << std::setprecision(FLOAT_PRECISION + 1) << acopy << " != " << bcopy; \
  } while (0)

inline Energy GetEnergy(const BackendModelPtr& m, const std::tuple<Primary, Secondary>& s) {
  return TotalEnergy(m, std::get<Primary>(s), std::get<Secondary>(s), nullptr).energy;
}

inline Energy GetEnergy(const BackendModelPtr& m, const std::string& r, const std::string& db) {
  return GetEnergy(m, {Primary::FromSeq(r), Secondary::FromDb(db)});
}

inline std::tuple<Energy, std::string> GetMfe(
    const BackendModelPtr& m, CtxCfg::MfeAlg alg, const Primary& r) {
  auto res = Ctx(m, CtxCfg{.mfe_alg = alg}).Fold(r, {});
  return {res.mfe.energy, mrna::BackendEnergyCfg(m).ToCtdString(res.tb.s, res.tb.ctd)};
}

inline std::tuple<Energy, std::string> GetMfe(
    const BackendModelPtr& m, CtxCfg::MfeAlg alg, const std::string& s) {
  return GetMfe(m, alg, Primary::FromSeq(s));
}

inline std::vector<subopt::SuboptResult> CheckSubopt(const BackendModelPtr& m,
    CtxCfg::SuboptAlg alg, const Primary& r, const std::vector<Energy>& energies) {
  const int n = static_cast<int>(energies.size());
  auto res = Ctx(m, CtxCfg{.subopt_alg = alg}).SuboptIntoVector(r, subopt::SuboptCfg{.strucs = n});
  for (int i = 0; i < n; ++i) EXPECT_EQ(res[i].energy, energies[i]);
  // for (int i = 0; i < n; ++i) fmt::print("E({}),\n", res[i].energy);
  return res;
}

inline std::vector<subopt::SuboptResult> CheckSubopt(const BackendModelPtr& m,
    CtxCfg::SuboptAlg alg, const std::string& s, const std::vector<Energy>& energies) {
  return CheckSubopt(m, alg, Primary::FromSeq(s), energies);
}

inline pfn::PfnResult GetPfn(const BackendModelPtr& m, CtxCfg::PfnAlg alg, const Primary& r) {
  return Ctx(m, CtxCfg{.pfn_alg = alg}).Pfn(r);
}

inline pfn::PfnResult GetPfn(const BackendModelPtr& m, CtxCfg::PfnAlg alg, const std::string& s) {
  return GetPfn(m, alg, Primary::FromSeq(s));
}

void CheckPfn(const PfnTables& got, const PfnTables& want);

std::tuple<Energy, Energy> GetPseudofree(
    const BackendModelPtr& m_base, const std::string& r, const std::string& db);

}  // namespace mrna

#endif  // TESTS_UTIL_H_
