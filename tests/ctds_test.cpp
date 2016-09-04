#include "constants.h"
#include "parsing.h"
#include "gtest/gtest.h"
#include "energy/energy_internal.h"
#include "energy/load_model.h"
#include "common_test.h"

namespace memerna {
namespace energy {

struct ctd_test_t {
  computed_t computed;
  internal::branch_ctd_t branch_ctds;
  std::deque<int> branches;
};

std::vector<ctd_test_t> CTD_TESTS;

void InitCtdsTest() {
  auto em = LoadEnergyModelFromDataDir(ENERGY_MODEL_PATH);
  secondary_t secondary = parsing::ParseDotBracketSecondary("GAAACAGAAAAUGGAAACCAGAAACA", "(...).((...).(...)).(...).");
  CTD_TESTS = std::vector<ctd_test_t>{
      {{}, {}, {}},
      {{parsing::ParseDotBracketSecondary("A", "."), {CTD_NA}, 0}, {}, {}},
      {{parsing::ParseDotBracketSecondary("AG", ".."), {CTD_NA, CTD_NA}, 0}, {}, {}},
      {{parsing::ParseDotBracketSecondary("GUA", "..."), {CTD_NA, CTD_NA, CTD_NA}, 0}, {}, {}},
      {{parsing::ParseDotBracketSecondary("GUAC", "...."), {CTD_NA, CTD_NA, CTD_NA, CTD_NA}, 0}, {}, {}},
      // 3' dangle inside the branch.
      {
          {
              parsing::ParseDotBracketSecondary("GAAAC", "(...)"),
              {CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_3_DANGLE}, 0
          },
          {{CTD_3_DANGLE, em.dangle3[G][A][C]}},
          {4}
      },
      {{secondary, std::vector<Ctd>(secondary.r.size(), CTD_NA), 0}, {}, {}},
      {
          {
              secondary,
              {
                  CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
                  CTD_RIGHT_MISMATCH_COAX_WITH_NEXT, CTD_NA,
                  CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
                  CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
                  CTD_NA, CTD_RIGHT_MISMATCH_COAX_WITH_PREV, CTD_NA, CTD_NA,
                  CTD_NA, CTD_NA, CTD_NA
              }, 0
          },
          {
              {CTD_UNUSED, 0}, {CTD_RIGHT_MISMATCH_COAX_WITH_NEXT, em.MismatchCoaxial(C, A, A, G)},
              {CTD_RIGHT_MISMATCH_COAX_WITH_PREV, em.MismatchCoaxial(C, A, A, G)}
          },
          {0, 6, 20}
      },
      {
          {
              secondary,
              {
                  CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
                  CTD_NA, CTD_FLUSH_COAX_WITH_PREV,
                  CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_5_DANGLE,
                  CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FLUSH_COAX_WITH_NEXT,
                  CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA
              }, 0
          },
          {
              {CTD_FLUSH_COAX_WITH_NEXT, em.stack[G][A][U][C]},
              {CTD_FLUSH_COAX_WITH_PREV, em.stack[G][A][U][C]},
              {CTD_5_DANGLE, em.dangle5[C][G][G]}
          },
          {18, 7, 13}
      },
      {
          {
              parsing::ParseDotBracketSecondary("GGAAACGAAACC", "((...)(...))"),
              {
                  CTD_NA, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FLUSH_COAX_WITH_NEXT,
                  CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FLUSH_COAX_WITH_PREV
              }, 0
          },
          {
              {CTD_UNUSED, 0}, {CTD_FLUSH_COAX_WITH_NEXT, em.stack[G][G][C][C]},
              {CTD_FLUSH_COAX_WITH_PREV, em.stack[G][G][C][C]}
          },
          {1, 6, 11}
      },
      {
          {
              parsing::ParseDotBracketSecondary("UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))"),
              {
                  CTD_NA, CTD_NA, CTD_NA, CTD_LEFT_MISMATCH_COAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
                  CTD_LEFT_MISMATCH_COAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
                  CTD_NA, CTD_FLUSH_COAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FLUSH_COAX_WITH_PREV
              }, 0
          },
          {
              {CTD_FLUSH_COAX_WITH_PREV, em.stack[U][C][G][A]},
              {CTD_LEFT_MISMATCH_COAX_WITH_NEXT, em.MismatchCoaxial(C, G, A, G)},
              {CTD_LEFT_MISMATCH_COAX_WITH_PREV, em.MismatchCoaxial(C, G, A, G)},
              {CTD_FLUSH_COAX_WITH_NEXT, em.stack[U][C][G][A]}
          },
          {24, 3, 9, 19}
      }
  };
}

class CtdsTest : public testing::TestWithParam<ctd_test_t> {
};

TEST_P(CtdsTest, BaseBranchBase) {
  // TODO unhack
  auto em = LoadEnergyModelFromDataDir(ENERGY_MODEL_PATH);
  auto ctd_test = GetParam();
  // Convert base representation to branch representation.
  internal::branch_ctd_t computed_branch_ctds;
  auto computed_energy = internal::GetBranchCtdsFromComputed(
      ctd_test.computed, em, ctd_test.branches, computed_branch_ctds);
  energy_t test_energy = 0;
  for (const auto& branch_ctd : ctd_test.branch_ctds) {
    // Make sure each branch energy is only represented once.
    if (branch_ctd.first == CTD_FLUSH_COAX_WITH_NEXT ||
        branch_ctd.first == CTD_LEFT_MISMATCH_COAX_WITH_NEXT ||
        branch_ctd.first == CTD_RIGHT_MISMATCH_COAX_WITH_NEXT)
      continue;
    test_energy += branch_ctd.second;
  }
  EXPECT_EQ(test_energy, computed_energy);
  EXPECT_EQ(ctd_test.branch_ctds, computed_branch_ctds);
  // Convert back again and make sure it's the same.
  std::vector<Ctd> previous_base_ctds = std::move(ctd_test.computed.base_ctds);
  ctd_test.computed.base_ctds.resize(previous_base_ctds.size(), CTD_NA);
  internal::AddBranchCtdsToComputed(ctd_test.computed, em, ctd_test.branches, computed_branch_ctds);
  EXPECT_EQ(previous_base_ctds, ctd_test.computed.base_ctds);
}

INSTANTIATE_TEST_CASE_P(CtdsTest, CtdsTest, testing::ValuesIn(CTD_TESTS));

}
}
