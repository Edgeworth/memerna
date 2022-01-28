// Copyright 2016 E.
#include <deque>
#include <functional>
#include <memory>
#include <tuple>
#include <utility>

#include "common_test.h"
#include "compute/energy/branch.h"
#include "compute/energy/model.h"
#include "gtest/gtest.h"
#include "model/base.h"
#include "model/ctd.h"
#include "model/model.h"
#include "model/primary.h"
#include "model/secondary.h"

namespace mrna::energy {

struct CtdTest {
  Primary r;
  Secondary s;
  Ctds ctd;
  BranchCtd branch_ctd;
  std::deque<int> branches;
};

std::function<CtdTest(const EnergyModel&)> CTD_TESTS[] = {[](const EnergyModel&) -> CtdTest {
                                                            return {{}, {}, {}, {}, {}};
                                                          },
    [](const EnergyModel&) -> CtdTest {
      return {Primary::FromString("A"), Secondary::FromDotBracket("."), {CTD_NA}, {}, {}};
    },
    [](const EnergyModel&) -> CtdTest {
      return {Primary::FromString("AG"), Secondary::FromDotBracket(".."), {CTD_NA, CTD_NA}, {}, {}};
    },
    [](const EnergyModel&) -> CtdTest {
      return {Primary::FromString("GUA"), Secondary::FromDotBracket("..."),
          {CTD_NA, CTD_NA, CTD_NA}, {}, {}};
    },
    [](const EnergyModel&) -> CtdTest {
      return {Primary::FromString("GUAC"), Secondary::FromDotBracket("...."),
          {CTD_NA, CTD_NA, CTD_NA, CTD_NA}, {}, {}};
    },
    // 3' dangle inside the branch.
    [](const EnergyModel& em) -> CtdTest {
      return {Primary::FromString("GAAAC"), Secondary::FromDotBracket("(...)"),
          {CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_3_DANGLE}, {{CTD_3_DANGLE, em.dangle3[G][A][C]}},
          {4}};
    },
    [](const EnergyModel&) -> CtdTest {
      return {Primary::FromString("GAAACAGAAAAUGGAAACCAGAAACA"),
          Secondary::FromDotBracket("(...).((...).(...)).(...)."), Ctds(26), {}, {}};
    },
    [](const EnergyModel& em) -> CtdTest {
      return {Primary::FromString("GAAACAGAAAAUGGAAACCAGAAACA"),
          Secondary::FromDotBracket("(...).((...).(...)).(...)."),
          {CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_RCOAX_WITH_NEXT, CTD_NA, CTD_NA,
              CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
              CTD_NA, CTD_RCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA},
          {{CTD_UNUSED, 0}, {CTD_RCOAX_WITH_NEXT, em.MismatchCoaxial(C, A, A, G)},
              {CTD_RCOAX_WITH_PREV, em.MismatchCoaxial(C, A, A, G)}},
          {0, 6, 20}};
    },
    [](const EnergyModel& em) -> CtdTest {
      return {Primary::FromString("GAAACAGAAAAUGGAAACCAGAAACA"),
          Secondary::FromDotBracket("(...).((...).(...)).(...)."),
          {CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_PREV, CTD_NA,
              CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_5_DANGLE, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
              CTD_FCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA},
          {{CTD_FCOAX_WITH_NEXT, em.stack[G][A][U][C]}, {CTD_FCOAX_WITH_PREV, em.stack[G][A][U][C]},
              {CTD_5_DANGLE, em.dangle5[C][G][G]}},
          {18, 7, 13}};
    },
    [](const EnergyModel& em) -> CtdTest {
      return {Primary::FromString("GGAAACGAAACC"), Secondary::FromDotBracket("((...)(...))"),
          {CTD_NA, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_NEXT, CTD_NA, CTD_NA,
              CTD_NA, CTD_NA, CTD_FCOAX_WITH_PREV},
          {{CTD_UNUSED, 0}, {CTD_FCOAX_WITH_NEXT, em.stack[G][G][C][C]},
              {CTD_FCOAX_WITH_PREV, em.stack[G][G][C][C]}},
          {1, 6, 11}};
    },
    [](const EnergyModel& em) -> CtdTest {
      return {Primary::FromString("UUAGAAACGCAAAGAGGUCCAAAGA"),
          Secondary::FromDotBracket("(..(...).(...).....(...))"),
          {CTD_NA, CTD_NA, CTD_NA, CTD_LCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
              CTD_LCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
              CTD_NA, CTD_FCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_PREV},
          {{CTD_FCOAX_WITH_PREV, em.stack[U][C][G][A]},
              {CTD_LCOAX_WITH_NEXT, em.MismatchCoaxial(C, G, A, G)},
              {CTD_LCOAX_WITH_PREV, em.MismatchCoaxial(C, G, A, G)},
              {CTD_FCOAX_WITH_NEXT, em.stack[U][C][G][A]}},
          {24, 3, 9, 19}};
    }};

class CtdsTest
    : public testing::TestWithParam<std::tuple<int, std::function<CtdTest(const EnergyModel&)>>> {};

TEST_P(CtdsTest, BaseBranchBase) {
  const auto& em = g_em[std::get<0>(GetParam())];
  auto ctd_test = std::get<1>(GetParam())(em);
  // Convert base representation to branch representation.
  BranchCtd computed_branch_ctd;
  auto computed_energy = AddBaseCtdsToBranchCtds(
      em, ctd_test.r, ctd_test.s, ctd_test.ctd, ctd_test.branches, &computed_branch_ctd);
  Energy test_energy = 0;
  for (const auto& branch_ctd : ctd_test.branch_ctd) {
    // Make sure each branch energy is only represented once.
    if (branch_ctd.first == CTD_FCOAX_WITH_NEXT || branch_ctd.first == CTD_LCOAX_WITH_NEXT ||
        branch_ctd.first == CTD_RCOAX_WITH_NEXT)
      continue;
    test_energy += branch_ctd.second;
  }
  EXPECT_EQ(test_energy, computed_energy);
  EXPECT_EQ(ctd_test.branch_ctd, computed_branch_ctd);
  // Convert back again and make sure it's the same.
  Ctds prev_ctd = std::move(ctd_test.ctd);
  ctd_test.ctd.reset(prev_ctd.size());
  AddBranchCtdsToBaseCtds(ctd_test.branches, computed_branch_ctd, &ctd_test.ctd);
  EXPECT_EQ(prev_ctd, ctd_test.ctd);
}

INSTANTIATE_TEST_SUITE_P(CtdsTest, CtdsTest,
    testing::Combine(testing::Range(0, NUM_TEST_MODELS), testing::ValuesIn(CTD_TESTS)));

}  // namespace mrna::energy
