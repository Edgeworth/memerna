// Copyright 2016 E.
#include "common_test.h"
#include "fold/brute_fold.h"
#include "gtest/gtest.h"
#include "model/parsing.h"

namespace mrna {
namespace fold {

TEST(BruteFold, GetBranchCounts) {
  EXPECT_EQ((std::vector<int>{2, 0}), internal::GetBranchCounts(DotBracketToPairs("()")));
  EXPECT_EQ((std::vector<int>{}), internal::GetBranchCounts(DotBracketToPairs("")));
  EXPECT_EQ((std::vector<int>{0}), internal::GetBranchCounts(DotBracketToPairs(".")));
  EXPECT_EQ(
      (std::vector<int>{2, 0, 2, 0, 2, 0}), internal::GetBranchCounts(DotBracketToPairs("()()()")));
  EXPECT_EQ((std::vector<int>{2, 1, 0, 1}), internal::GetBranchCounts(DotBracketToPairs("(())")));
  EXPECT_EQ((std::vector<int>{2, 3, 0, 3, 0, 3, 0, 3}),
      internal::GetBranchCounts(DotBracketToPairs("(()()())")));
  EXPECT_EQ((std::vector<int>{2, 1, 2, 0, 2, 0, 2, 1}),
      internal::GetBranchCounts(DotBracketToPairs("((()()))")));
}

}  // namespace fold
}  // namespace mrna
