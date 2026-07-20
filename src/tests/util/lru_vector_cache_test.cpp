// Copyright 2026 Eliot Courtney.
#include "util/lru_vector_cache.h"

#include <vector>

#include "gtest/gtest.h"

namespace mrna {

TEST(LruVectorCacheTest, EvictsLeastRecentlyUsedValues) {
  LruVectorCache<int> cache(/*max_size=*/1000, /*max_total=*/2);
  cache.Insert(10, std::vector<int>{1});
  cache.Insert(20, std::vector<int>{2});
  EXPECT_EQ(cache.Get(10), (std::vector<int>{1}));

  cache.Insert(30, std::vector<int>{3});
  EXPECT_TRUE(cache.Get(20).empty());
  EXPECT_EQ(cache.Get(10), (std::vector<int>{1}));
  EXPECT_EQ(cache.Get(30), (std::vector<int>{3}));
}

}  // namespace mrna
