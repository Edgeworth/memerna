// Copyright 2016 Eliot Courtney.
#include "model/base.h"

#include <memory>

#include "gtest/gtest.h"

namespace mrna {

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

}  // namespace mrna
