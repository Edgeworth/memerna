// Copyright 2016 Eliot Courtney.
#include "model/base.h"

#include <gtest/gtest.h>

#include <string>

namespace mrna {

TEST(BaseTest, CanPair) {
  EXPECT_TRUE(IsPair(A, U));
  EXPECT_TRUE(IsPair(U, A));
  EXPECT_TRUE(IsPair(G, U));
  EXPECT_TRUE(IsPair(U, G));
  EXPECT_TRUE(IsPair(G, C));
  EXPECT_TRUE(IsPair(C, G));
  EXPECT_FALSE(IsPair(A, A));
  EXPECT_FALSE(IsPair(C, C));
  EXPECT_FALSE(IsPair(G, G));
  EXPECT_FALSE(IsPair(U, U));
  EXPECT_FALSE(IsPair(A, C));
  EXPECT_FALSE(IsPair(C, A));
  EXPECT_FALSE(IsPair(A, G));
  EXPECT_FALSE(IsPair(G, A));
  EXPECT_FALSE(IsPair(C, U));
  EXPECT_FALSE(IsPair(U, C));
}

TEST(BaseTest, IsAuGu) {
  EXPECT_TRUE(IsAuGuPair(A, U));
  EXPECT_TRUE(IsAuGuPair(U, A));
  EXPECT_TRUE(IsAuGuPair(G, U));
  EXPECT_TRUE(IsAuGuPair(U, G));
  EXPECT_FALSE(IsAuGuPair(G, C));
  EXPECT_FALSE(IsAuGuPair(C, G));
  EXPECT_FALSE(IsAuGuPair(A, A));
  EXPECT_FALSE(IsAuGuPair(C, C));
  EXPECT_FALSE(IsAuGuPair(G, G));
  EXPECT_FALSE(IsAuGuPair(U, U));
  EXPECT_FALSE(IsAuGuPair(A, C));
  EXPECT_FALSE(IsAuGuPair(C, A));
  EXPECT_FALSE(IsAuGuPair(A, G));
  EXPECT_FALSE(IsAuGuPair(G, A));
  EXPECT_FALSE(IsAuGuPair(C, U));
  EXPECT_FALSE(IsAuGuPair(U, C));
}

}  // namespace mrna
