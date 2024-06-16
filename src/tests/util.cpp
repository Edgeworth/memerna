// Copyright 2024 Eliot Courtney.
#include "tests/util.h"

#include "api/ctx/backend.h"
namespace mrna {

void CheckPfn(const PfnTables& got, const PfnTables& want) {
  const int N = static_cast<int>(want.p.size());
  EXPECT_EQ(want.prob.size(), want.p.size());
  EXPECT_EQ(got.p.size(), want.p.size());
  EXPECT_EQ(got.prob.size(), want.p.size());
  EXPECT_REL_EQ(got.q, want.q);

  for (int st = 0; st < N; ++st) {
    for (int en = 0; en < N; ++en) {
      EXPECT_REL_EQ(got.p[st][en], want.p[st][en]);
      EXPECT_REL_EQ(got.prob[st][en], want.prob[st][en]);
    }
  }
}

std::tuple<Energy, Energy> GetPseudofree(
    const BackendModelPtr& m_base, const std::string& r, const std::string& db) {
  const auto paired_mul = E(100.0);
  const auto unpaired_mul = E(10.0);
  std::vector<Energy> pf_paired(r.size(), E(0.0));
  std::vector<Energy> pf_unpaired(r.size(), E(0.0));
  Energy extra_from_pseudofree = ZERO_E;
  for (int i = 0; i < int(r.size()); ++i) {
    pf_paired[i] = paired_mul * (i + 1);
    pf_unpaired[i] = unpaired_mul * (i + 1);
    if (db[i] == '.') {
      extra_from_pseudofree += pf_unpaired[i];
    } else {
      extra_from_pseudofree += pf_paired[i];
    }
  }

  auto m_pf = CloneBackend(m_base);
  LoadPseudofreeEnergy(m_pf, pf_paired, pf_unpaired);
  auto energy = GetEnergy(m_pf, {Primary::FromSeq(r), Secondary::FromDb(db)});
  return {energy, extra_from_pseudofree};
}

}  // namespace mrna
