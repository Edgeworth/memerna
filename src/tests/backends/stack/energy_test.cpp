// Copyright 2016 Eliot Courtney.
#include "model/energy.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>

#include "gtest/gtest.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "tests/init.h"

namespace mrna::md::stack {

class EnergyTestStack : public testing::TestWithParam<int> {
 public:
  static Energy GetEnergy(const std::string& r, const std::string& db) {
    return GetEnergy({Primary::FromSeq(r), Secondary::FromDb(db)});
  }

  static Energy GetEnergy(const std::tuple<Primary, Secondary>& s) {
    return stack_ms[GetParam()]
        ->TotalEnergy(std::get<Primary>(s), std::get<Secondary>(s), nullptr)
        .energy;
  }
};

TEST_P(EnergyTestStack, MultiloopEnergy) {
  const auto& m = stack_ms[GetParam()];
  EXPECT_EQ(m->multiloop_a + 4 * m->multiloop_b, m->MultiloopInitiation(4));
}

TEST_P(EnergyTestStack, NNDBHairpinLoopExamples) {
  const auto& m = stack_ms[GetParam()];

  EXPECT_EQ(m->stack[C][A][U][G] + m->stack[A][C][G][U] + m->stack[C][A][U][G] + m->au_penalty +
          m->terminal[A][A][A][U] + m->HairpinInitiation(6) + m->penultimate_stack[C][A][U][G] +
          m->penultimate_stack[U][G][C][A],
      GetEnergy(kNNDBHairpin1));
  EXPECT_EQ(m->stack[C][A][U][G] + m->stack[A][C][G][U] + m->stack[C][A][U][G] + m->au_penalty +
          m->terminal[A][G][G][U] + m->hairpin_gg_first_mismatch + m->HairpinInitiation(5) +
          m->penultimate_stack[C][A][U][G] + m->penultimate_stack[U][G][C][A],
      GetEnergy(kNNDBHairpin2));

  if (m->hairpin.contains("CCGAGG")) {
    EXPECT_EQ(m->stack[C][A][U][G] + m->stack[A][C][G][U] + m->stack[C][C][G][G] +
            m->hairpin["CCGAGG"] + m->penultimate_stack[C][C][G][G] +
            m->penultimate_stack[U][G][C][A],
        GetEnergy(kNNDBHairpin3));
  }

  EXPECT_EQ(m->stack[C][A][U][G] + m->stack[A][C][G][U] + m->stack[C][A][U][G] + m->au_penalty +
          m->terminal[A][C][C][U] + m->HairpinInitiation(6) + m->hairpin_all_c_a * 6 +
          m->hairpin_all_c_b + m->penultimate_stack[C][A][U][G] + m->penultimate_stack[U][G][C][A],
      GetEnergy(kNNDBHairpin4));
  EXPECT_EQ(m->stack[C][G][C][G] + m->stack[G][G][C][C] + m->stack[G][G][U][C] + m->gu_penalty +
          m->terminal[G][G][G][U] + m->hairpin_gg_first_mismatch + m->HairpinInitiation(5) +
          m->hairpin_special_gu_closure + m->penultimate_stack[G][G][U][C] +
          m->penultimate_stack[C][G][C][G],
      GetEnergy(kNNDBHairpin5));
}

TEST_P(EnergyTestStack, NNDBBulgeLoopExamples) {
  const auto& m = stack_ms[GetParam()];

  EXPECT_EQ(m->stack[G][C][G][C] + m->stack[C][C][G][G] + m->BulgeInitiation(1) +
          m->bulge_special_c + m->stack[C][G][C][G] + m->HairpinInitiation(3) - E(R * T * log(3)) +
          m->penultimate_stack[C][G][C][G] + m->penultimate_stack[G][C][G][C],
      GetEnergy(kNNDBBulge1));

  EXPECT_EQ(m->stack[G][A][U][C] + m->au_penalty + m->BulgeInitiation(3) + m->HairpinInitiation(3) +
          m->penultimate_stack[U][C][G][A] + m->penultimate_stack[G][A][U][C],
      GetEnergy(kNNDBBulge2));
}

TEST_P(EnergyTestStack, NNDBInternalLoopExamples) {
  const auto& m = stack_ms[GetParam()];

  EXPECT_EQ(m->stack[C][A][U][G] + m->stack[C][G][C][G] + m->InternalLoopInitiation(5) +
          std::min(m->internal_asym, NINIO_MAX_ASYM) + m->internal_2x3_mismatch[A][G][G][U] +
          m->internal_2x3_mismatch[G][G][A][C] + m->au_penalty + m->internal_au_penalty +
          m->HairpinInitiation(3) + m->penultimate_stack[C][G][C][G] +
          m->penultimate_stack[C][G][C][G] + m->penultimate_stack[C][A][U][G] +
          m->penultimate_stack[U][G][C][A],
      GetEnergy(kNNDBInternal2x3));
  EXPECT_EQ(m->stack[C][A][U][G] + m->stack[C][G][C][G] + m->internal_2x2[A][G][A][C][G][G][A][U] +
          m->au_penalty + m->HairpinInitiation(3) + m->penultimate_stack[C][G][C][G] +
          m->penultimate_stack[C][G][C][G] + m->penultimate_stack[C][A][U][G] +
          m->penultimate_stack[U][G][C][A],
      GetEnergy(kNNDBInternal2x2));
  EXPECT_EQ(m->stack[C][A][U][G] + m->stack[C][G][C][G] + m->InternalLoopInitiation(6) +
          std::min(4 * m->internal_asym, NINIO_MAX_ASYM) + m->au_penalty + m->internal_au_penalty +
          m->HairpinInitiation(3) + m->penultimate_stack[C][G][C][G] +
          m->penultimate_stack[C][G][C][G] + m->penultimate_stack[C][A][U][G] +
          m->penultimate_stack[U][G][C][A],
      GetEnergy(kNNDBInternal1x5));
}

TEST_P(EnergyTestStack, BaseCases) {
  const auto& m = stack_ms[GetParam()];

  EXPECT_EQ(m->au_penalty + m->stack[G][A][U][C] + m->hairpin_init[3] +
          m->penultimate_stack[G][A][U][C] + m->penultimate_stack[U][C][G][A],
      GetEnergy("GAAAAUC", "((...))"));
  EXPECT_EQ(m->au_penalty + m->gu_penalty + m->stack[G][A][U][U] + m->hairpin_init[3] +
          m->penultimate_stack[G][A][U][U] + m->penultimate_stack[U][U][G][A],
      GetEnergy("GAAAAUU", "((...))"));
  EXPECT_EQ(m->au_penalty * 2 + m->HairpinInitiation(3) +
          std::min({ZERO_E, m->terminal[U][A][A][A], m->dangle3[U][A][A], m->dangle5[U][A][A]}),
      GetEnergy("AAAAAUA", ".(...)."));
  EXPECT_EQ(m->au_penalty * 2 + m->HairpinInitiation(3), GetEnergy("AAAAU", "(...)"));
  EXPECT_EQ(m->stack[G][C][G][C] + m->stack[C][U][A][G] + m->BulgeInitiation(1) +
          m->stack[U][G][C][A] + m->HairpinInitiation(3) + m->penultimate_stack[U][G][C][A] +
          m->penultimate_stack[G][C][G][C],
      GetEnergy(kBulge1));
  EXPECT_EQ(m->InternalLoopInitiation(5) + std::min(m->internal_asym, NINIO_MAX_ASYM) +
          m->internal_au_penalty + m->au_penalty * 2 + m->internal_2x3_mismatch[A][G][A][U] +
          m->internal_2x3_mismatch[C][A][A][G] + m->HairpinInitiation(3),
      GetEnergy(kInternal1));
}

INSTANTIATE_TEST_SUITE_P(EnergyModelTests, EnergyTestStack, testing::Range(0, NUM_TEST_MODELS));

#if ENERGY_PRECISION == 1

TEST(EnergyTestStack, T04) {
  auto m = stack_t04;

  EXPECT_EQ(E(8.8), m->HairpinInitiation(87));
  EXPECT_EQ(E(6.8), m->BulgeInitiation(57));
  EXPECT_EQ(E(4.6), m->InternalLoopInitiation(67));
}

#elif ENERGY_PRECISION == 2

TEST(EnergyTestStack, T04) {
  auto m = stack_t04;

  EXPECT_EQ(E(8.85), m->HairpinInitiation(87));
  EXPECT_EQ(E(6.79), m->BulgeInitiation(57));
  EXPECT_EQ(E(4.57), m->InternalLoopInitiation(67));
}

#endif

}  // namespace mrna::md::stack
