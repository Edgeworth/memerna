// Copyright 2016 Eliot Courtney.
#include "util/splaymap.h"

#include <cstdint>
#include <limits>
#include <random>
#include <set>

#include "gtest/gtest.h"
#include "util/error.h"

namespace mrna {

TEST(SplayMapTest, SmallCase1) {
  SplayMap<int, int> h;
  for (int i = 0; i < 10; ++i) {
    EXPECT_EQ(0u, h.Size());
    EXPECT_FALSE(h.Delete(1));
    EXPECT_FALSE(h.Find(1));
    EXPECT_TRUE(h.Insert(1, 1));
    EXPECT_EQ(1u, h.Size());
    EXPECT_TRUE(h.Find(1));
    EXPECT_EQ(1, h.Get());
    EXPECT_FALSE(h.Find(2));
    EXPECT_TRUE(h.Insert(2, 2));
    EXPECT_TRUE(h.Find(2));
    EXPECT_EQ(2, h.Get());
    EXPECT_TRUE(h.Find(1));
    EXPECT_EQ(1, h.Get());
    EXPECT_EQ(2u, h.Size());
    EXPECT_TRUE(h.Delete(2));
    EXPECT_EQ(1u, h.Size());
    EXPECT_FALSE(h.Delete(2));
    EXPECT_EQ(1u, h.Size());
    EXPECT_TRUE(h.Delete(1));
    EXPECT_EQ(0u, h.Size());
    EXPECT_FALSE(h.Delete(0));
  }
}

TEST(SplayMapTest, SmallCase2) {
  SplayMap<int, int> h;
  EXPECT_EQ(0u, h.Size());
  EXPECT_TRUE(h.Insert(1, 1));
  EXPECT_TRUE(h.Insert(2, 2));
  EXPECT_EQ(2u, h.Size());
  EXPECT_EQ(2, h.Get());
  EXPECT_FALSE(h.Insert(1, 3));
  EXPECT_EQ(1, h.Get());
  EXPECT_EQ(2u, h.Size());
}

class SplayMapRandomTest : public testing::TestWithParam<int_fast32_t> {};

TEST_P(SplayMapRandomTest, CompareAgainstMap) {
  const int NUM_TRIES = 1000;
  std::mt19937 eng(GetParam());
  SplayMap<int, int> h;
  std::set<int> s;
  std::vector<int> keys;
  auto val_dist = std::uniform_int_distribution<int>(
      std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
  auto case_dist = std::uniform_int_distribution<int>(0, 4);
  for (int i = 0; i < NUM_TRIES; ++i) {
    EXPECT_EQ(s.size(), h.Size());
    auto h_keys = h.Keys();
    EXPECT_TRUE(std::equal(s.begin(), s.end(), h_keys.begin(), h_keys.end()));
    int key = 0;
    switch (case_dist(eng)) {
    case 0:
      // Check existing
      if (!keys.empty()) {
        key = keys[std::uniform_int_distribution<std::size_t>(0, keys.size() - 1)(eng)];
        EXPECT_TRUE(s.count(key));
        EXPECT_TRUE(h.Find(key));
        EXPECT_EQ(key, h.Get());
      }
      break;
    case 1:
      // Check non-existing
      // NOLINTNEXTLINE
      do { key = val_dist(eng); } while (s.count(key));
      EXPECT_FALSE(h.Find(key));
      break;
    case 2: {
      // Insert case
      int val = val_dist(eng);
      if (s.count(val))
        EXPECT_FALSE(h.Insert(val, val));
      else
        EXPECT_TRUE(h.Insert(val, val));
      s.insert(val);
      break;
    }
    case 3:
      // Delete existing
      if (!keys.empty()) {
        auto idx = std::uniform_int_distribution<std::size_t>(0, keys.size() - 1)(eng);
        std::swap(keys[idx], keys.back());
        key = keys.back();
        keys.pop_back();
        EXPECT_TRUE(s.count(key));
        EXPECT_TRUE(h.Delete(key));
      }
      break;
    case 4:
      // Delete non-existing
      // NOLINTNEXTLINE
      do { key = val_dist(eng); } while (s.count(key));
      EXPECT_FALSE(h.Delete(key));
      break;
    default: bug();
    }
  }
}

INSTANTIATE_TEST_SUITE_P(
    SplayMapRandomTest, SplayMapRandomTest, testing::Range(int_fast32_t(0), int_fast32_t(5)));

}  // namespace mrna
