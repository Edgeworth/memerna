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
#include "base.h"

#include "gtest/gtest.h"

namespace memerna {

TEST(BaseTest, CanPair) {
  EXPECT_TRUE(CanPair(A, U));
  EXPECT_TRUE(CanPair(U, A));
  EXPECT_TRUE(CanPair(G, U));
  EXPECT_TRUE(CanPair(U, G));
  EXPECT_TRUE(CanPair(G, C));
  EXPECT_TRUE(CanPair(C, G));
  EXPECT_FALSE(CanPair(A, A));
  EXPECT_FALSE(CanPair(C, C));
  EXPECT_FALSE(CanPair(G, G));
  EXPECT_FALSE(CanPair(U, U));
  EXPECT_FALSE(CanPair(A, C));
  EXPECT_FALSE(CanPair(C, A));
  EXPECT_FALSE(CanPair(A, G));
  EXPECT_FALSE(CanPair(G, A));
  EXPECT_FALSE(CanPair(C, U));
  EXPECT_FALSE(CanPair(U, C));
}

TEST(BaseTest, IsAuGu) {
  EXPECT_TRUE(IsAuGu(A, U));
  EXPECT_TRUE(IsAuGu(U, A));
  EXPECT_TRUE(IsAuGu(G, U));
  EXPECT_TRUE(IsAuGu(U, G));
  EXPECT_FALSE(IsAuGu(G, C));
  EXPECT_FALSE(IsAuGu(C, G));
  EXPECT_FALSE(IsAuGu(A, A));
  EXPECT_FALSE(IsAuGu(C, C));
  EXPECT_FALSE(IsAuGu(G, G));
  EXPECT_FALSE(IsAuGu(U, U));
  EXPECT_FALSE(IsAuGu(A, C));
  EXPECT_FALSE(IsAuGu(C, A));
  EXPECT_FALSE(IsAuGu(A, G));
  EXPECT_FALSE(IsAuGu(G, A));
  EXPECT_FALSE(IsAuGu(C, U));
  EXPECT_FALSE(IsAuGu(U, C));
}
}  // namespace memerna
