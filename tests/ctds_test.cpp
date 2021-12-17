// Copyright 2016 Eliot Courtney.
#include "common_test.h"
#include "compute/energy/internal.h"
#include "compute/energy/load_model.h"
#include "gtest/gtest.h"
#include "model/parsing.h"

namespace mrna {
namespace energy {

struct ctd_test_t {
  computed_t computed;
  internal::branch_ctd_t branch_ctds;
  std::deque<int> branches;
};

std::function<ctd_test_t(const EnergyModel&)> CTD_TESTS[] = {[](const EnergyModel&) -> ctd_test_t {
                                                               return {{}, {}, {}};
                                                             },
    [](const EnergyModel&) -> ctd_test_t {
      return {{ParseDotBracketSecondary("A", "."), {CTD_NA}, 0}, {}, {}};
    },
    [](const EnergyModel&) -> ctd_test_t {
      return {{ParseDotBracketSecondary("AG", ".."), {CTD_NA, CTD_NA}, 0}, {}, {}};
    },
    [](const EnergyModel&) -> ctd_test_t {
      return {{ParseDotBracketSecondary("GUA", "..."), {CTD_NA, CTD_NA, CTD_NA}, 0}, {}, {}};
    },
    [](const EnergyModel&) -> ctd_test_t {
      return {
          {ParseDotBracketSecondary("GUAC", "...."), {CTD_NA, CTD_NA, CTD_NA, CTD_NA}, 0}, {}, {}};
    },
    // 3' dangle inside the branch.
    [](const EnergyModel& em) -> ctd_test_t {
      return {{ParseDotBracketSecondary("GAAAC", "(...)"),
                  {CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_3_DANGLE}, 0},
          {{CTD_3_DANGLE, em.dangle3[G][A][C]}}, {4}};
    },
    [](const EnergyModel&) -> ctd_test_t {
      return {{ParseDotBracketSecondary("GAAACAGAAAAUGGAAACCAGAAACA", "(...).((...).(...)).(...)."),
                  std::vector<Ctd>(26, CTD_NA), 0},
          {}, {}};
    },
    [](const EnergyModel& em) -> ctd_test_t {
      return {
          {ParseDotBracketSecondary("GAAACAGAAAAUGGAAACCAGAAACA", "(...).((...).(...)).(...)."),
              {CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_RCOAX_WITH_NEXT, CTD_NA,
                  CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
                  CTD_NA, CTD_NA, CTD_RCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA},
              0},
          {{CTD_UNUSED, 0}, {CTD_RCOAX_WITH_NEXT, em.MismatchCoaxial(C, A, A, G)},
              {CTD_RCOAX_WITH_PREV, em.MismatchCoaxial(C, A, A, G)}},
          {0, 6, 20}};
    },
    [](const EnergyModel& em) -> ctd_test_t {
      return {
          {ParseDotBracketSecondary("GAAACAGAAAAUGGAAACCAGAAACA", "(...).((...).(...)).(...)."),
              {CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_PREV, CTD_NA,
                  CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_5_DANGLE, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
                  CTD_FCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA},
              0},
          {{CTD_FCOAX_WITH_NEXT, em.stack[G][A][U][C]}, {CTD_FCOAX_WITH_PREV, em.stack[G][A][U][C]},
              {CTD_5_DANGLE, em.dangle5[C][G][G]}},
          {18, 7, 13}};
    },
    [](const EnergyModel& em) -> ctd_test_t {
      return {{ParseDotBracketSecondary("GGAAACGAAACC", "((...)(...))"),
                  {CTD_NA, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_NEXT, CTD_NA,
                      CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_PREV},
                  0},
          {{CTD_UNUSED, 0}, {CTD_FCOAX_WITH_NEXT, em.stack[G][G][C][C]},
              {CTD_FCOAX_WITH_PREV, em.stack[G][G][C][C]}},
          {1, 6, 11}};
    },
    [](const EnergyModel& em) -> ctd_test_t {
      return {{ParseDotBracketSecondary("UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))"),
                  {CTD_NA, CTD_NA, CTD_NA, CTD_LCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
                      CTD_NA, CTD_LCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
                      CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
                      CTD_FCOAX_WITH_PREV},
                  0},
          {{CTD_FCOAX_WITH_PREV, em.stack[U][C][G][A]},
              {CTD_LCOAX_WITH_NEXT, em.MismatchCoaxial(C, G, A, G)},
              {CTD_LCOAX_WITH_PREV, em.MismatchCoaxial(C, G, A, G)},
              {CTD_FCOAX_WITH_NEXT, em.stack[U][C][G][A]}},
          {24, 3, 9, 19}};
    }};

class CtdsTest : public testing::TestWithParam<
                     std::tuple<EnergyModelPtr, std::function<ctd_test_t(const EnergyModel&)>>> {};

TEST_P(CtdsTest, BaseBranchBase) {
  const auto& em = *std::get<0>(GetParam());
  auto ctd_test = std::get<1>(GetParam())(em);
  // Convert base representation to branch representation.
  internal::branch_ctd_t computed_branch_ctds;
  auto computed_energy = internal::GetBranchCtdsFromComputed(
      ctd_test.computed, em, ctd_test.branches, computed_branch_ctds);
  energy_t test_energy = 0;
  for (const auto& branch_ctd : ctd_test.branch_ctds) {
    // Make sure each branch energy is only represented once.
    if (branch_ctd.first == CTD_FCOAX_WITH_NEXT || branch_ctd.first == CTD_LCOAX_WITH_NEXT ||
        branch_ctd.first == CTD_RCOAX_WITH_NEXT)
      continue;
    test_energy += branch_ctd.second;
  }
  EXPECT_EQ(test_energy, computed_energy);
  EXPECT_EQ(ctd_test.branch_ctds, computed_branch_ctds);
  // Convert back again and make sure it's the same.
  std::vector<Ctd> previous_base_ctds = std::move(ctd_test.computed.base_ctds);
  ctd_test.computed.base_ctds.resize(previous_base_ctds.size(), CTD_NA);
  internal::AddBranchCtdsToComputed(ctd_test.computed, ctd_test.branches, computed_branch_ctds);
  EXPECT_EQ(previous_base_ctds, ctd_test.computed.base_ctds);
}

INSTANTIATE_TEST_SUITE_P(CtdsTest, CtdsTest,
    testing::Combine(testing::Values(LoadRandomEnergyModel(0), LoadRandomEnergyModel(1),
                         LoadRandomEnergyModel(2), LoadRandomEnergyModel(3)),
        testing::ValuesIn(CTD_TESTS)));

}  // namespace energy
}  // namespace mrna
