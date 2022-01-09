// Copyright 2016 E.
#include "compute/brute.h"

#include "common_test.h"
#include "gtest/gtest.h"

namespace mrna {

TEST(BruteFold, GetBranchCounts) {
  EXPECT_EQ((std::vector<int>{2, 0}), internal::GetBranchCounts(DotBracketToSecondary("()")));
  EXPECT_EQ((std::vector<int>{}), internal::GetBranchCounts(DotBracketToSecondary("")));
  EXPECT_EQ((std::vector<int>{0}), internal::GetBranchCounts(DotBracketToSecondary(".")));
  EXPECT_EQ((std::vector<int>{2, 0, 2, 0, 2, 0}),
      internal::GetBranchCounts(DotBracketToSecondary("()()()")));
  EXPECT_EQ(
      (std::vector<int>{2, 1, 0, 1}), internal::GetBranchCounts(DotBracketToSecondary("(())")));
  EXPECT_EQ((std::vector<int>{2, 3, 0, 3, 0, 3, 0, 3}),
      internal::GetBranchCounts(DotBracketToSecondary("(()()())")));
  EXPECT_EQ((std::vector<int>{2, 1, 2, 0, 2, 0, 2, 1}),
      internal::GetBranchCounts(DotBracketToSecondary("((()()))")));
}

}  // namespace mrna
