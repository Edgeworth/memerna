// Copyright 2024 Eliot Courtney.
#include "api/subopt/subopt.h"

#include <string>
#include <tuple>
#include <vector>

#include "api/ctx/ctx_cfg.h"
#include "gtest/gtest.h"
#include "model/energy.h"
#include "model/primary.h"
#include "tests/init.h"
#include "tests/util.h"

namespace mrna {

class SuboptTestT12 : public testing::TestWithParam<std::tuple<int, CtxCfg::SuboptAlg>> {
 public:
  static std::vector<subopt::SuboptResult> Subopt(
      const BackendModelPtr& m, const std::string& s, const std::vector<Energy>& energies) {
    return CheckSubopt(m, std::get<1>(GetParam()), s, energies);
  }

  static std::vector<subopt::SuboptResult> Subopt(
      const BackendModelPtr& m, const Primary& r, const std::vector<Energy>& energies) {
    return CheckSubopt(m, std::get<1>(GetParam()), r, energies);
  }
};

#if ENERGY_PRECISION == 1

GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(SuboptTestT12);

#elif ENERGY_PRECISION == 2

TEST_P(SuboptTestT12, T12P2) {
  auto [i, alg] = GetParam();
  const auto& m = t12_ms[i];
  if (!Contains(CtxCfg::SuboptAlgsForBackend(m), alg)) return;

  Subopt(m, "CCUCCGGG",
      {
          E(-0.56),
          E(0.00),
          E(0.84),
          E(1.20),
          E(1.40),
          E(2.14),
          E(2.40),
          E(2.70),
          E(3.33),
          E(3.40),
      });
  Subopt(m, "CGGAAACGG",
      {
          E(-0.56),
          E(0.00),
          E(0.74),
          E(1.40),
          E(1.60),
          E(2.60),
          E(2.80),
          E(3.10),
          E(3.90),
      });
  Subopt(m, "UGCAAAGCAA",
      {
          E(-0.43),
          E(0.00),
          E(0.28),
          E(0.37),
          E(0.48),
          E(1.98),
          E(1.98),
          E(2.19),
          E(2.90),
          E(2.99),
      });
  Subopt(m, "GGGGAAACCCC",
      {
          E(-4.38),
          E(-4.22),
          E(-3.42),
          E(-2.42),
          E(-2.42),
          E(-2.32),
          E(-2.12),
          E(-1.92),
          E(-1.37),
          E(-1.37),
      });
  Subopt(m, "GGGGAAACCCC",
      {
          E(-4.38),
          E(-4.22),
          E(-3.42),
          E(-2.42),
          E(-2.42),
          E(-2.32),
          E(-2.12),
          E(-1.92),
          E(-1.37),
          E(-1.37),
      });
  Subopt(m, "CUUAUAGUUAAGG",
      {
          E(-1.24),
          E(0.00),
          E(0.06),
          E(0.79),
          E(1.64),
          E(1.84),
          E(1.92),
          E(2.09),
          E(2.54),
          E(2.57),
      });
  Subopt(m, "CCGAAGGGGCUGCGGCG",
      {
          E(-3.94),
          E(-2.38),
          E(-2.32),
          E(-2.22),
          E(-2.18),
          E(-1.50),
          E(-1.48),
          E(-1.42),
          E(-1.36),
          E(-1.12),
      });
  Subopt(m, "GCCAAGGCCCCACCCGGA",
      {
          E(-2.88),
          E(-2.04),
          E(-1.18),
          E(-0.92),
          E(-0.92),
          E(-0.78),
          E(-0.72),
          E(-0.66),
          E(-0.49),
          E(-0.38),
      });
  Subopt(m, "GGCCGAUGGCAGCGAUAGC",
      {
          E(-4.90),
          E(-3.10),
          E(-2.98),
          E(-2.80),
          E(-2.68),
          E(-2.34),
          E(-1.60),
          E(-1.60),
          E(-1.40),
          E(-1.40),
      });
  Subopt(m, "CUGAAACUGGAAACAGAAAUG",
      {
          E(-2.19),
          E(-2.19),
          E(-1.59),
          E(-1.09),
          E(0.00),
          E(0.11),
          E(0.31),
          E(0.41),
          E(0.44),
          E(0.49),
      });

  // Too slow for brute force:
  if (alg == CtxCfg::SuboptAlg::BRUTE) return;
  Subopt(m, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA",
      {
          E(-5.19),
          E(-5.07),
          E(-5.07),
          E(-4.85),
          E(-4.79),
          E(-4.63),
          E(-4.61),
          E(-4.61),
          E(-4.59),
          E(-4.58),
      });
  Subopt(m, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU",
      {
          E(-12.94),
          E(-12.89),
          E(-12.83),
          E(-12.69),
          E(-12.53),
          E(-12.49),
          E(-12.38),
          E(-12.35),
          E(-12.33),
          E(-12.30),
      });
  Subopt(m, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC",
      {
          E(-5.59),
          E(-5.58),
          E(-5.39),
          E(-5.38),
          E(-5.38),
          E(-5.18),
          E(-5.09),
          E(-4.89),
          E(-4.85),
          E(-4.84),
      });
  Subopt(m, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG",
      {
          E(-12.03),
          E(-11.73),
          E(-11.53),
          E(-11.23),
          E(-11.13),
          E(-10.93),
          E(-10.93),
          E(-10.83),
          E(-10.63),
          E(-10.63),
      });
  Subopt(m, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG",
      {
          E(-7.28),
          E(-7.18),
          E(-6.32),
          E(-6.32),
          E(-5.74),
          E(-5.70),
          E(-5.70),
          E(-5.64),
          E(-5.60),
          E(-5.60),
      });
  Subopt(m, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG",
      {
          E(-2.83),
          E(-2.47),
          E(-2.43),
          E(-2.23),
          E(-2.23),
          E(-2.05),
          E(-2.03),
          E(-1.94),
          E(-1.94),
          E(-1.92),
      });
  Subopt(m, "CCGGGCCAGCCCGCUCCUACGGGGGGUC",
      {
          E(-11.83),
          E(-11.00),
          E(-10.69),
          E(-10.69),
          E(-10.14),
          E(-10.05),
          E(-10.01),
          E(-9.93),
          E(-9.71),
          E(-9.61),
      });
  Subopt(m, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG",
      {
          E(-7.49),
          E(-7.19),
          E(-7.19),
          E(-6.89),
          E(-6.49),
          E(-6.29),
          E(-6.19),
          E(-5.99),
          E(-5.92),
          E(-5.89),
      });
  Subopt(m, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC",
      {
          E(-3.18),
          E(-2.98),
          E(-2.98),
          E(-2.87),
          E(-2.78),
          E(-2.68),
          E(-2.68),
          E(-2.48),
          E(-2.48),
          E(-2.45),
      });
  Subopt(m, "CGCAGGGUCGGACCCGGGAGAACCGCGA",
      {
          E(-6.36),
          E(-6.04),
          E(-5.65),
          E(-5.42),
          E(-5.32),
          E(-5.30),
          E(-5.28),
          E(-5.28),
          E(-5.26),
          E(-5.12),
      });
  Subopt(m, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA",
      {
          E(-5.93),
          E(-5.73),
          E(-5.63),
          E(-5.59),
          E(-5.54),
          E(-5.43),
          E(-5.41),
          E(-5.39),
          E(-5.37),
          E(-5.35),
      });
  Subopt(m, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC",
      {
          E(-12.05),
          E(-12.05),
          E(-12.05),
          E(-12.05),
          E(-11.97),
          E(-11.97),
          E(-11.97),
          E(-11.97),
          E(-11.97),
          E(-11.97),
      });
  Subopt(m, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC",
      {
          E(-3.95),
          E(-3.75),
          E(-3.70),
          E(-3.60),
          E(-3.50),
          E(-3.40),
          E(-3.28),
          E(-3.15),
          E(-3.06),
          E(-3.06),
      });
  Subopt(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC",
      {
          E(-5.93),
          E(-5.93),
          E(-5.87),
          E(-5.53),
          E(-5.53),
          E(-5.47),
          E(-5.44),
          E(-5.23),
          E(-5.23),
          E(-5.17),
      });
  Subopt(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
      {
          E(-28.29),
          E(-27.91),
          E(-27.91),
          E(-27.74),
          E(-27.71),
          E(-27.71),
          E(-27.64),
          E(-27.64),
          E(-27.62),
          E(-27.58),
      });
  Subopt(m, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU",
      {
          E(-5.39),
          E(-5.06),
          E(-4.96),
          E(-4.16),
          E(-3.94),
          E(-3.86),
          E(-3.76),
          E(-3.73),
          E(-3.68),
          E(-3.66),
      });
  Subopt(m, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG",
      {
          E(-15.06),
          E(-14.86),
          E(-14.20),
          E(-14.00),
          E(-13.70),
          E(-13.56),
          E(-13.53),
          E(-13.53),
          E(-13.50),
          E(-13.49),
      });
  Subopt(m, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC",
      {
          E(-4.18),
          E(-4.00),
          E(-3.91),
          E(-3.90),
          E(-3.90),
          E(-3.87),
          E(-3.84),
          E(-3.74),
          E(-3.71),
          E(-3.70),
      });
  Subopt(m, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG",
      {
          E(-3.25),
          E(-3.05),
          E(-3.05),
          E(-3.05),
          E(-2.76),
          E(-2.64),
          E(-2.55),
          E(-2.48),
          E(-2.47),
          E(-2.44),
      });
  Subopt(m, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG",
      {
          E(-2.29),
          E(-1.98),
          E(-1.98),
          E(-1.83),
          E(-1.82),
          E(-1.79),
          E(-1.78),
          E(-1.78),
          E(-1.71),
          E(-1.71),
      });
  Subopt(m, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG",
      {
          E(-8.53),
          E(-8.53),
          E(-8.23),
          E(-8.23),
          E(-8.14),
          E(-8.03),
          E(-8.03),
          E(-7.98),
          E(-7.89),
          E(-7.68),
      });
  Subopt(m, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU",
      {
          E(-19.99),
          E(-19.99),
          E(-19.99),
          E(-19.64),
          E(-19.64),
          E(-19.64),
          E(-19.59),
          E(-19.49),
          E(-19.49),
          E(-19.44),
      });
  Subopt(m, std::get<Primary>(k16sHSapiens3),
      {
          E(-88.99),
          E(-88.99),
          E(-88.96),
          E(-88.96),
          E(-88.89),
          E(-88.89),
          E(-88.86),
          E(-88.86),
          E(-88.79),
          E(-88.79),
      });
}

#endif

INSTANTIATE_TEST_SUITE_P(SuboptTest, SuboptTestT12,
    testing::Combine(
        testing::Range(0, NUM_T12_MODELS), testing::ValuesIn(EnumValues<CtxCfg::SuboptAlg>())));

}  // namespace mrna
