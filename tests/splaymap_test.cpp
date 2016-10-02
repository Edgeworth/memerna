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
#include "splaymap.h"
#include "gtest/gtest.h"

namespace memerna {

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

class SplayMapRandomTest : public testing::TestWithParam<int_fast32_t> {
};

TEST_P(SplayMapRandomTest, CompareAgainstMap) {
  const int NUM_TRIES = 10000;
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
        do key = val_dist(eng);
        while (s.count(key));
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
        do key = val_dist(eng);
        while (s.count(key));
        EXPECT_FALSE(h.Delete(key));
        break;
      default:
        verify_expr(false, "bug");
    }
  }
}

INSTANTIATE_TEST_CASE_P(SplayMapRandomTest, SplayMapRandomTest,
    testing::Range(int_fast32_t(0), int_fast32_t(10)));

}
