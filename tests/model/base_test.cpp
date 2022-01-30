// Copyright 2016 E.
#include "model/base.h"

#include <memory>

#include "gtest/gtest.h"

namespace mrna {

TEST(BaseTest, CanPair) {
  EXPECT_TRUE(BasePair(A, U));
  EXPECT_TRUE(BasePair(U, A));
  EXPECT_TRUE(BasePair(G, U));
  EXPECT_TRUE(BasePair(U, G));
  EXPECT_TRUE(BasePair(G, C));
  EXPECT_TRUE(BasePair(C, G));
  EXPECT_FALSE(BasePair(A, A));
  EXPECT_FALSE(BasePair(C, C));
  EXPECT_FALSE(BasePair(G, G));
  EXPECT_FALSE(BasePair(U, U));
  EXPECT_FALSE(BasePair(A, C));
  EXPECT_FALSE(BasePair(C, A));
  EXPECT_FALSE(BasePair(A, G));
  EXPECT_FALSE(BasePair(G, A));
  EXPECT_FALSE(BasePair(C, U));
  EXPECT_FALSE(BasePair(U, C));
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

}  // namespace mrna
