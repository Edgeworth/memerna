// Copyright 2016 Eliot Courtney.
#include "model/base.h"

#include "gtest/gtest.h"

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
  EXPECT_TRUE(IsAuPair(A, U));
  EXPECT_TRUE(IsAuPair(U, A));

  EXPECT_TRUE(IsGuPair(G, U));
  EXPECT_TRUE(IsGuPair(U, G));

  EXPECT_FALSE(IsAuPair(G, C));
  EXPECT_FALSE(IsGuPair(G, C));

  EXPECT_FALSE(IsAuPair(C, G));
  EXPECT_FALSE(IsGuPair(C, G));

  EXPECT_FALSE(IsAuPair(A, A));
  EXPECT_FALSE(IsGuPair(A, A));

  EXPECT_FALSE(IsAuPair(C, C));
  EXPECT_FALSE(IsGuPair(C, C));

  EXPECT_FALSE(IsAuPair(G, G));
  EXPECT_FALSE(IsGuPair(G, G));

  EXPECT_FALSE(IsAuPair(U, U));
  EXPECT_FALSE(IsGuPair(U, U));

  EXPECT_FALSE(IsAuPair(A, C));
  EXPECT_FALSE(IsGuPair(A, C));

  EXPECT_FALSE(IsAuPair(C, A));
  EXPECT_FALSE(IsGuPair(C, A));

  EXPECT_FALSE(IsAuPair(A, G));
  EXPECT_FALSE(IsGuPair(A, G));

  EXPECT_FALSE(IsAuPair(G, A));
  EXPECT_FALSE(IsGuPair(G, A));

  EXPECT_FALSE(IsAuPair(C, U));
  EXPECT_FALSE(IsGuPair(C, U));

  EXPECT_FALSE(IsAuPair(U, C));
  EXPECT_FALSE(IsGuPair(U, C));
}

}  // namespace mrna
