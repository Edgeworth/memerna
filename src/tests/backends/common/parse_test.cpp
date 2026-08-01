// Copyright 2026 Eliot Courtney.
#include "backends/common/parse.h"

#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "gtest/gtest.h"
#include "model/base.h"

namespace mrna::md {

TEST(ParseTest, FixedMapRejectsDuplicateEntries) {
  Energy output[4][4][4][4]{};
  EXPECT_THROW(Parse4Map(std::istringstream("AAAA 0\nAAAA 1\n"), output), std::runtime_error);
}

TEST(ParseTest, FixedMapClearsAndAcceptsSparseEntries) {
  Energy output[4][4][4][4];
  for (auto& a : output)
    for (auto& b : a)
      for (auto& c : b)
        for (auto& energy : c) energy = E(1.0);
  EXPECT_NO_THROW(Parse4Map(std::istringstream("AGGU 2\n"), output));
  EXPECT_EQ(E(2.0), output[A][G][G][U]);
  EXPECT_EQ(ZERO_E, output[A][A][A][U]);
}

TEST(ParseTest, VectorClearsAndAcceptsSparseEntries) {
  Energy output[3]{E(1.0), E(1.0), E(1.0)};
  EXPECT_NO_THROW(ParseVec(std::istringstream("1 2\n"), output));
  EXPECT_EQ(E(2.0), output[1]);
  EXPECT_EQ(ZERO_E, output[2]);
}

TEST(ParseTest, NMapRejectsDuplicateEntries) {
  std::unordered_map<std::string, Energy> output;
  EXPECT_THROW(ParseNMap(std::istringstream("CAACG 0\nCAACG 1\n"), &output), std::runtime_error);
}

TEST(ParseTest, NMapAcceptsCustomEntries) {
  std::unordered_map<std::string, Energy> output;
  EXPECT_NO_THROW(ParseNMap(std::istringstream("custom 2\n"), &output));
  EXPECT_EQ(E(2.0), output.at("custom"));
}

}  // namespace mrna::md
