// Copyright 2016 Eliot Courtney.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <deque>
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

#include "common_test.h"
#include "compute/energy/branch.h"
#include "compute/energy/energy.h"
#include "compute/energy/t04/branch.h"
#include "compute/energy/t04/model.h"
#include "compute/energy/t04/precomp.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"

namespace mrna::erg {

class T22ModelTest : public testing::TestWithParam<int> {
 public:
  std::tuple<Primary, Secondary> kNNDBHairpin1 = ParseSeqDb("CACAAAAAAAUGUG", "((((......))))");
  std::tuple<Primary, Secondary> kNNDBHairpin2 = ParseSeqDb("CACAGGAAGUGUG", "((((.....))))");
  std::tuple<Primary, Secondary> kNNDBHairpin3 = ParseSeqDb("CACCCGAGGGUG", "((((....))))");
  std::tuple<Primary, Secondary> kNNDBHairpin4 = ParseSeqDb("CACACCCCCCUGUG", "((((......))))");
  std::tuple<Primary, Secondary> kNNDBHairpin5 = ParseSeqDb("CGGGGGAAGUCCG", "((((.....))))");
  std::tuple<Primary, Secondary> kNNDBBulge1 = ParseSeqDb("GCCCGAAACGGC", "(((.(...))))");
  std::tuple<Primary, Secondary> kNNDBBulge2 = ParseSeqDb("GAACAGAAACUC", "((...(...)))");
  std::tuple<Primary, Secondary> kNNDBInternal2x3 =
      ParseSeqDb("CAGACGAAACGGAGUG", "((..((...))...))");
  std::tuple<Primary, Secondary> kNNDBInternal1x5 =
      ParseSeqDb("CAGCGAAACGGAAAGUG", "((.((...)).....))");
  std::tuple<Primary, Secondary> kNNDBInternal2x2 =
      ParseSeqDb("CAGACGAAACGGAUG", "((..((...))..))");

  std::tuple<Primary, Secondary> kBulge1 = ParseSeqDb("GCUCGAAACAGC", "(((.(...))))");
  std::tuple<Primary, Secondary> kInternal1 = ParseSeqDb("AGAGAAACAAAU", "(..(...)...)");

  static Energy GetEnergy(const std::string& r, const std::string& db) { return E(0); }

  static Energy GetEnergy(const std::tuple<Primary, Secondary>& s) { return E(0); }
};

TEST_P(T22ModelTest, MultiloopEnergy) {}

TEST_P(T22ModelTest, NNDBHairpinLoopExamples) {}

TEST_P(T22ModelTest, NNDBBulgeLoopExamples) {}

TEST_P(T22ModelTest, NNDBInternalLoopExamples) {}

TEST_P(T22ModelTest, BaseCases) {}

INSTANTIATE_TEST_SUITE_P(EnergyModelTests, T22ModelTest, testing::Range(0, NUM_TEST_MODELS));

}  // namespace mrna::erg
