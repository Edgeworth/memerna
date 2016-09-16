#include "minmaxheap.h"
#include "gtest/gtest.h"

namespace memerna {

TEST(MinMaxHeapTest, ZeroCase) {
  MinMaxHeap<int> h;
  EXPECT_EQ(0u, h.Size());
  h.Insert(1);
  EXPECT_EQ(1u, h.Size());
  EXPECT_EQ(1, h.Min());
  EXPECT_EQ(1, h.Max());
  h.PopMin();
  EXPECT_EQ(0u, h.Size());
  h.Insert(1);
  EXPECT_EQ(1u, h.Size());
  EXPECT_EQ(1, h.Min());
  EXPECT_EQ(1, h.Max());
  h.PopMax();
  EXPECT_EQ(0u, h.Size());
}

TEST(MinMaxHeapTest, SmallCase) {
  MinMaxHeap<int> h;
  EXPECT_EQ(0u, h.Size());
  h.Insert(1);
  EXPECT_EQ(1, h.Min());
  EXPECT_EQ(1, h.Max());
  h.Insert(2);
  EXPECT_EQ(2u, h.Size());
  EXPECT_EQ(1, h.Min());
  EXPECT_EQ(2, h.Max());
  h.PopMax();
  EXPECT_EQ(1, h.Min());
  EXPECT_EQ(1, h.Max());
  h.PopMin();
  EXPECT_EQ(0u, h.Size());
}

class MinMaxRandomTest : public testing::TestWithParam<int_fast32_t> {
};

TEST_P(MinMaxRandomTest, CompareAgainstSet) {
  std::mt19937 eng(GetParam());
  MinMaxHeap<int> h;
  std::multiset<int> s;
  int len = std::uniform_int_distribution<int>(5, 100000)(eng);
  auto valdist = std::uniform_int_distribution<int>(-50, 50);
  for (int i = 0; i < len; ++i) {
    int val = valdist(eng);
    h.Insert(val);
    s.insert(val);

    EXPECT_EQ(*s.begin(), h.Min());
    EXPECT_EQ(*(--s.end()), h.Max());
    EXPECT_EQ(s.size(), h.Size());
  }
  //h.PrintTreeDebug();
  for (int i = 0; i < len; ++i) {
    EXPECT_EQ(*s.begin(), h.Min());
    EXPECT_EQ(*(--s.end()), h.Max());
    EXPECT_EQ(s.size(), h.Size());

    if (eng() % 2) {
      //printf("Popping min\n");
      h.PopMin();
      s.erase(s.begin());
    } else {
      //printf("Popping max\n");
      h.PopMax();
      s.erase(--s.end());
    }
    //h.PrintTreeDebug();
  }
  EXPECT_EQ(0u, h.Size());
}

INSTANTIATE_TEST_CASE_P(MinMaxRandomTest, MinMaxRandomTest,
    testing::Range(int_fast32_t(0), int_fast32_t(10))); // TODO chagne back

}
