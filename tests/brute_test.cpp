// Copyright 2016, E.
//
// This file is part of memerna.
//
// memerna is free software: you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// memerna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with memerna.
// If not, see <http://www.gnu.org/licenses/>.
#include "gtest/gtest.h"
#include "common_test.h"
#include "fold/brute_fold.h"
#include "parsing.h"

namespace memerna {
namespace fold {

TEST(BruteFold, GetBranchCounts) {
  EXPECT_EQ((std::vector<int>{2, 0}), internal::GetBranchCounts(parsing::DotBracketToPairs("()")));
  EXPECT_EQ((std::vector<int>{}), internal::GetBranchCounts(parsing::DotBracketToPairs("")));
  EXPECT_EQ((std::vector<int>{0}), internal::GetBranchCounts(parsing::DotBracketToPairs(".")));
  EXPECT_EQ((std::vector<int>{2, 0, 2, 0, 2, 0}),
      internal::GetBranchCounts(parsing::DotBracketToPairs("()()()")));
  EXPECT_EQ((std::vector<int>{2, 1, 0, 1}),
      internal::GetBranchCounts(parsing::DotBracketToPairs("(())")));
  EXPECT_EQ((std::vector<int>{2, 3, 0, 3, 0, 3, 0, 3}),
      internal::GetBranchCounts(parsing::DotBracketToPairs("(()()())")));
  EXPECT_EQ((std::vector<int>{2, 1, 2, 0, 2, 0, 2, 1}),
      internal::GetBranchCounts(parsing::DotBracketToPairs("((()()))")));
}
}
}
