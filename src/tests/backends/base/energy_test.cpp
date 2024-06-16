// Copyright 2016 Eliot Courtney.
#include "model/energy.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <deque>
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

#include "backends/base/energy/precomp.h"
#include "backends/common/base/branch.h"
#include "backends/common/base/model_base.h"
#include "backends/common/branch.h"
#include "gtest/gtest.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "tests/init.h"

namespace mrna::md::base {

class EnergyTestBase : public testing::TestWithParam<int> {
 public:
  static Energy GetEnergy(const std::string& r, const std::string& db) {
    return GetEnergy({Primary::FromSeq(r), Secondary::FromDb(db)});
  }

  static Energy GetEnergy(const std::tuple<Primary, Secondary>& s) {
    return base_ms[GetParam()]
        ->TotalEnergy(std::get<Primary>(s), std::get<Secondary>(s), nullptr)
        .energy;
  }
};

TEST_P(EnergyTestBase, MultiloopEnergy) {
  auto m = base_ms[GetParam()];
  EXPECT_EQ(m->multiloop_hack_a + 4 * m->multiloop_hack_b, m->MultiloopInitiation(4));
}

TEST_P(EnergyTestBase, NNDBHairpinLoopExamples) {
  auto m = base_ms[GetParam()];

  EXPECT_EQ(m->stack[C][A][U][G] + m->stack[A][C][G][U] + m->stack[C][A][U][G] + m->au_penalty +
          m->terminal[A][A][A][U] + m->HairpinInitiation(6),
      GetEnergy(kNNDBHairpin1));
  EXPECT_EQ(m->stack[C][A][U][G] + m->stack[A][C][G][U] + m->stack[C][A][U][G] + m->au_penalty +
          m->terminal[A][G][G][U] + m->hairpin_gg_first_mismatch + m->HairpinInitiation(5),
      GetEnergy(kNNDBHairpin2));

  if (m->hairpin.contains("CCGAGG")) {
    EXPECT_EQ(
        m->stack[C][A][U][G] + m->stack[A][C][G][U] + m->stack[C][C][G][G] + m->hairpin["CCGAGG"],
        GetEnergy(kNNDBHairpin3));
  }

  EXPECT_EQ(m->stack[C][A][U][G] + m->stack[A][C][G][U] + m->stack[C][A][U][G] + m->au_penalty +
          m->terminal[A][C][C][U] + m->HairpinInitiation(6) + m->hairpin_all_c_a * 6 +
          m->hairpin_all_c_b,
      GetEnergy(kNNDBHairpin4));
  EXPECT_EQ(m->stack[C][G][C][G] + m->stack[G][G][C][C] + m->stack[G][G][U][C] + m->gu_penalty +
          m->terminal[G][G][G][U] + m->hairpin_gg_first_mismatch + m->HairpinInitiation(5) +
          m->hairpin_special_gu_closure,
      GetEnergy(kNNDBHairpin5));

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin1)), m);
    EXPECT_EQ(m->au_penalty + m->terminal[A][A][A][U] + m->HairpinInitiation(6), pc.Hairpin(3, 10));
  }

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin2)), m);
    EXPECT_EQ(m->au_penalty + m->terminal[A][G][G][U] + m->hairpin_gg_first_mismatch +
            m->HairpinInitiation(5),
        pc.Hairpin(3, 9));
  }

  if (m->hairpin.contains("CCGAGG")) {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin3)), m);
    EXPECT_EQ(m->hairpin["CCGAGG"], pc.Hairpin(3, 8));
  }

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin4)), m);
    EXPECT_EQ(m->au_penalty + m->terminal[A][C][C][U] + m->HairpinInitiation(6) +
            m->hairpin_all_c_a * 6 + m->hairpin_all_c_b,
        pc.Hairpin(3, 10));
  }

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin5)), m);
    EXPECT_EQ(m->gu_penalty + m->terminal[G][G][G][U] + m->hairpin_gg_first_mismatch +
            m->HairpinInitiation(5) + m->hairpin_special_gu_closure,
        pc.Hairpin(3, 9));
  }
}

TEST_P(EnergyTestBase, NNDBBulgeLoopExamples) {
  auto m = base_ms[GetParam()];

  EXPECT_EQ(m->stack[G][C][G][C] + m->stack[C][C][G][G] + m->BulgeInitiation(1) +
          m->bulge_special_c + m->stack[C][G][C][G] + m->HairpinInitiation(3) - E(R * T * log(3)),
      GetEnergy(kNNDBBulge1));

  EXPECT_EQ(m->stack[G][A][U][C] + m->au_penalty + m->BulgeInitiation(3) + m->HairpinInitiation(3),
      GetEnergy(kNNDBBulge2));
}

TEST_P(EnergyTestBase, NNDBInternalLoopExamples) {
  auto m = base_ms[GetParam()];

  EXPECT_EQ(m->stack[C][A][U][G] + m->stack[C][G][C][G] + m->InternalLoopInitiation(5) +
          std::min(m->internal_asym, NINIO_MAX_ASYM) + m->internal_2x3_mismatch[A][G][G][U] +
          m->internal_2x3_mismatch[G][G][A][C] + m->au_penalty + m->internal_au_penalty +
          m->HairpinInitiation(3),
      GetEnergy(kNNDBInternal2x3));
  EXPECT_EQ(m->stack[C][A][U][G] + m->stack[C][G][C][G] + m->internal_2x2[A][G][A][C][G][G][A][U] +
          m->au_penalty + m->HairpinInitiation(3),
      GetEnergy(kNNDBInternal2x2));
  EXPECT_EQ(m->stack[C][A][U][G] + m->stack[C][G][C][G] + m->InternalLoopInitiation(6) +
          std::min(4 * m->internal_asym, NINIO_MAX_ASYM) + m->au_penalty + m->internal_au_penalty +
          m->HairpinInitiation(3),
      GetEnergy(kNNDBInternal1x5));
}

TEST_P(EnergyTestBase, BaseCases) {
  auto m = base_ms[GetParam()];

  EXPECT_EQ(
      m->au_penalty + m->stack[G][A][U][C] + m->hairpin_init[3], GetEnergy("GAAAAUC", "((...))"));
  EXPECT_EQ(m->au_penalty + m->gu_penalty + m->stack[G][A][U][U] + m->hairpin_init[3],
      GetEnergy("GAAAAUU", "((...))"));
  EXPECT_EQ(m->au_penalty * 2 + m->HairpinInitiation(3) +
          std::min(ZERO_E,
              std::min(
                  m->terminal[U][A][A][A], std::min(m->dangle3[U][A][A], m->dangle5[U][A][A]))),
      GetEnergy("AAAAAUA", ".(...)."));
  EXPECT_EQ(m->au_penalty * 2 + m->HairpinInitiation(3), GetEnergy("AAAAU", "(...)"));
  EXPECT_EQ(m->stack[G][C][G][C] + m->stack[C][U][A][G] + m->BulgeInitiation(1) +
          m->stack[U][G][C][A] + m->HairpinInitiation(3),
      GetEnergy(kBulge1));
  EXPECT_EQ(m->InternalLoopInitiation(5) + std::min(m->internal_asym, NINIO_MAX_ASYM) +
          m->internal_au_penalty + m->au_penalty * 2 + m->internal_2x3_mismatch[A][G][A][U] +
          m->internal_2x3_mismatch[C][A][A][G] + m->HairpinInitiation(3),
      GetEnergy(kInternal1));
}

INSTANTIATE_TEST_SUITE_P(EnergyModelTests, EnergyTestBase, testing::Range(0, NUM_TEST_MODELS));

#if ENERGY_PRECISION == 1

TEST(EnergyTestBase, T04) {
  auto m = base_t04;

  EXPECT_EQ(E(8.8), m->HairpinInitiation(87));
  EXPECT_EQ(E(6.8), m->BulgeInitiation(57));
  EXPECT_EQ(E(4.6), m->InternalLoopInitiation(67));

  const Precomp pc(Primary::FromSeq("GGGGAAACCCC"), m);
  EXPECT_EQ(E(-2.1 - 0.4 - 1.6), pc.min_mismatch_coax);
  EXPECT_EQ(E(-3.4), pc.min_flush_coax);
  EXPECT_EQ(E(-2.6), pc.min_twoloop_not_stack);

  Energy augubranch[4][4] = {{E(-0.6), E(-0.6), E(-0.6), E(0.5 - 0.6)},
      {E(-0.6), E(-0.6), E(-0.6), E(-0.6)}, {E(-0.6), E(-0.6), E(-0.6), E(0.5 - 0.6)},
      {E(0.5 - 0.6), E(-0.6), E(0.5 - 0.6), E(-0.6)}};
  EXPECT_EQ(sizeof(augubranch), sizeof(pc.augubranch));
  EXPECT_EQ(0, std::memcmp(augubranch, pc.augubranch, sizeof(augubranch)));
}

#elif ENERGY_PRECISION == 2

TEST(EnergyTestBase, T04) {
  auto m = base_t04;

  EXPECT_EQ(E(8.85), m->HairpinInitiation(87));
  EXPECT_EQ(E(6.79), m->BulgeInitiation(57));
  EXPECT_EQ(E(4.57), m->InternalLoopInitiation(67));

  const Precomp pc(Primary::FromSeq("GGGGAAACCCC"), m);
  EXPECT_EQ(E(-2.10 - 0.40 - 1.60), pc.min_mismatch_coax);
  EXPECT_EQ(E(-3.42), pc.min_flush_coax);
  EXPECT_EQ(E(-2.60), pc.min_twoloop_not_stack);

  Energy augubranch[4][4] = {{E(-0.60), E(-0.60), E(-0.60), E(0.50 - 0.60)},
      {E(-0.60), E(-0.60), E(-0.60), E(-0.60)}, {E(-0.60), E(-0.60), E(-0.60), E(0.50 - 0.60)},
      {E(0.50 - 0.60), E(-0.60), E(0.50 - 0.60), E(-0.60)}};
  EXPECT_EQ(sizeof(augubranch), sizeof(pc.augubranch));
  EXPECT_EQ(0, std::memcmp(augubranch, pc.augubranch, sizeof(augubranch)));
}

#endif

struct CtdTest {
  Primary r;
  Secondary s;
  Ctds ctd;
  BranchCtd branch_ctd;
  std::deque<int> branches;
};

std::function<CtdTest(const Model::Ptr&)> CTD_TESTS[] = {
    [](const Model::Ptr&) -> CtdTest { return {{}, {}, {}, {}, {}}; },
    [](const Model::Ptr&) -> CtdTest {
      return {Primary::FromSeq("A"), Secondary::FromDb("."), Ctds{CTD_NA}, {}, {}};
    },
    [](const Model::Ptr&) -> CtdTest {
      return {Primary::FromSeq("AG"), Secondary::FromDb(".."), Ctds{CTD_NA, CTD_NA}, {}, {}};
    },
    [](const Model::Ptr&) -> CtdTest {
      return {
          Primary::FromSeq("GUA"), Secondary::FromDb("..."), Ctds{CTD_NA, CTD_NA, CTD_NA}, {}, {}};
    },
    [](const Model::Ptr&) -> CtdTest {
      return {Primary::FromSeq("GUAC"), Secondary::FromDb("...."),
          Ctds{CTD_NA, CTD_NA, CTD_NA, CTD_NA}, {}, {}};
    },
    // 3' dangle inside the branch.
    [](const Model::Ptr& m) -> CtdTest {
      return {Primary::FromSeq("GAAAC"), Secondary::FromDb("(...)"),
          Ctds{CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_3_DANGLE}, {{CTD_3_DANGLE, m->dangle3[G][A][C]}},
          {4}};
    },
    [](const Model::Ptr&) -> CtdTest {
      return {Primary::FromSeq("GAAACAGAAAAUGGAAACCAGAAACA"),
          Secondary::FromDb("(...).((...).(...)).(...)."), Ctds(26), {}, {}};
    },
    [](const Model::Ptr& m) -> CtdTest {
      return {Primary::FromSeq("GAAACAGAAAAUGGAAACCAGAAACA"),
          Secondary::FromDb("(...).((...).(...)).(...)."),
          Ctds{CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_RC_WITH_NEXT, CTD_NA, CTD_NA,
              CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
              CTD_NA, CTD_RC_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA},
          {{CTD_UNUSED, ZERO_E}, {CTD_RC_WITH_NEXT, m->MismatchCoaxial(C, A, A, G)},
              {CTD_RC_WITH_PREV, m->MismatchCoaxial(C, A, A, G)}},
          {0, 6, 20}};
    },
    [](const Model::Ptr& m) -> CtdTest {
      return {Primary::FromSeq("GAAACAGAAAAUGGAAACCAGAAACA"),
          Secondary::FromDb("(...).((...).(...)).(...)."),
          Ctds{CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_PREV, CTD_NA,
              CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_5_DANGLE, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
              CTD_FCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA},
          {{CTD_FCOAX_WITH_NEXT, m->stack[G][A][U][C]}, {CTD_FCOAX_WITH_PREV, m->stack[G][A][U][C]},
              {CTD_5_DANGLE, m->dangle5[C][G][G]}},
          {18, 7, 13}};
    },
    [](const Model::Ptr& m) -> CtdTest {
      return {Primary::FromSeq("GGAAACGAAACC"), Secondary::FromDb("((...)(...))"),
          Ctds{CTD_NA, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_NEXT, CTD_NA,
              CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_PREV},
          {{CTD_UNUSED, ZERO_E}, {CTD_FCOAX_WITH_NEXT, m->stack[G][G][C][C]},
              {CTD_FCOAX_WITH_PREV, m->stack[G][G][C][C]}},
          {1, 6, 11}};
    },
    [](const Model::Ptr& m) -> CtdTest {
      return {Primary::FromSeq("UUAGAAACGCAAAGAGGUCCAAAGA"),
          Secondary::FromDb("(..(...).(...).....(...))"),
          Ctds{CTD_NA, CTD_NA, CTD_NA, CTD_LCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
              CTD_LCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
              CTD_NA, CTD_FCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_PREV},
          {{CTD_FCOAX_WITH_PREV, m->stack[U][C][G][A]},
              {CTD_LCOAX_WITH_NEXT, m->MismatchCoaxial(C, G, A, G)},
              {CTD_LCOAX_WITH_PREV, m->MismatchCoaxial(C, G, A, G)},
              {CTD_FCOAX_WITH_NEXT, m->stack[U][C][G][A]}},
          {24, 3, 9, 19}};
    }};

class CtdsTestBase
    : public testing::TestWithParam<std::tuple<int, std::function<CtdTest(const Model::Ptr&)>>> {};

TEST_P(CtdsTestBase, BaseBranchBase) {
  auto m = base_ms[std::get<0>(GetParam())];
  auto ctd_test = std::get<1>(GetParam())(m);
  // Convert base representation to branch representation.
  BranchCtd computed_branch_ctd;
  auto computed_energy = AddBaseCtdsToBranchCtds(
      *m, ctd_test.r, ctd_test.s, ctd_test.ctd, ctd_test.branches, &computed_branch_ctd);
  Energy test_energy = ZERO_E;
  for (const auto& branch_ctd : ctd_test.branch_ctd) {
    // Make sure each branch energy is only represented once.
    if (branch_ctd.first == CTD_FCOAX_WITH_NEXT || branch_ctd.first == CTD_LCOAX_WITH_NEXT ||
        branch_ctd.first == CTD_RC_WITH_NEXT)
      continue;
    test_energy += branch_ctd.second;
  }
  EXPECT_EQ(test_energy, computed_energy);
  EXPECT_EQ(ctd_test.branch_ctd, computed_branch_ctd);
  // Convert back again and make sure it's the same.
  const Ctds prev_ctd = std::move(ctd_test.ctd);
  ctd_test.ctd.reset(prev_ctd.size());
  AddBranchCtdsToBaseCtds(ctd_test.branches, computed_branch_ctd, &ctd_test.ctd);
  EXPECT_EQ(prev_ctd, ctd_test.ctd);
}

INSTANTIATE_TEST_SUITE_P(CtdsTest, CtdsTestBase,
    testing::Combine(testing::Range(0, NUM_TEST_MODELS), testing::ValuesIn(CTD_TESTS)));

}  // namespace mrna::md::base
