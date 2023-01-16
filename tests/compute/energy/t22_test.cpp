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
#include "compute/energy/common/branch.h"
#include "compute/energy/common/t04like/branch.h"
#include "compute/energy/energy.h"
#include "compute/energy/t04/model.h"
#include "compute/energy/t04/precomp.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"

namespace mrna::erg {

Energy GetEnergy(const t22::Model::Ptr& em, const std::string& r, const std::string& db) {
  return em->TotalEnergy(Primary::FromSeq(r), Secondary::FromDb(db), nullptr).energy;
}

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
};

TEST_P(T22ModelTest, MultiloopEnergy) {}

TEST_P(T22ModelTest, NNDBHairpinLoopExamples) {}

TEST_P(T22ModelTest, NNDBBulgeLoopExamples) {}

TEST_P(T22ModelTest, NNDBInternalLoopExamples) {}

TEST_P(T22ModelTest, BaseCases) {}

#if ENERGY_PRECISION == 2

TEST(T22P2ModelTest, T22Tests) {
  // TODO(0): Uncomment.
  // auto em = t22p2;

  // EXPECT_EQ(E(8.85), em->HairpinInitiation(87));
  // EXPECT_EQ(E(6.79), em->BulgeInitiation(57));
  // EXPECT_EQ(E(4.57), em->InternalLoopInitiation(67));

  // // Example from https://doi.org/10.1093/nar/gkac261
  // EXPECT_EQ(E(0), GetEnergy(em, "UGUCGAUACCCUGUCGAUA", "((((((((...))))))))"));
  // EXPECT_EQ(E(0), GetEnergy(em, "UAGGUCAGCCCCUGGUCUA", "((((((((...))))))))"));
}

#endif

INSTANTIATE_TEST_SUITE_P(EnergyModelTests, T22ModelTest, testing::Range(0, NUM_TEST_MODELS));

}  // namespace mrna::erg
