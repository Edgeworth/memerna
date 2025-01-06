// Copyright 2016 Eliot Courtney.
#include "model/energy.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <deque>
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

#include "backends/base/energy/precomp.h"
#include "backends/common/base/branch.h"
#include "backends/common/base/model_base.h"
#include "gtest/gtest.h"
#include "model/base.h"
#include "model/branch.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "tests/init.h"
#include "tests/util.h"

namespace mrna::md::base {

class EnergyTestBase : public testing::TestWithParam<int> {
 public:
  static Energy GetEnergy(const std::string& r, const std::string& db) {
    return GetEnergy({Primary::FromSeq(r), Secondary::FromDb(db)});
  }

  static Energy GetEnergy(const std::tuple<Primary, Secondary>& s) {
    return base_ms[GetParam()]
        ->TotalEnergy(std::get<Primary>(s), std::get<Secondary>(s), nullptr)
        .energy;
  }
};

TEST_P(EnergyTestBase, MultiloopEnergy) {
  const auto& m = base_ms[GetParam()];
  EXPECT_EQ(m->multiloop_hack_a + 4 * m->multiloop_hack_b, m->MultiloopInitiation(4));
}

TEST_P(EnergyTestBase, NNDBHairpinLoopExamples) {
  const auto& m = base_ms[GetParam()];

  EXPECT_EQ(m->stack[C][A][U][G] + m->stack[A][C][G][U] + m->stack[C][A][U][G] + m->au_penalty +
          m->terminal[A][A][A][U] + m->HairpinInitiation(6),
      GetEnergy(kNNDBHairpin1));
  EXPECT_EQ(m->stack[C][A][U][G] + m->stack[A][C][G][U] + m->stack[C][A][U][G] + m->au_penalty +
          m->terminal[A][G][G][U] + m->hairpin_gg_first_mismatch + m->HairpinInitiation(5),
      GetEnergy(kNNDBHairpin2));

  if (m->hairpin.contains("CCGAGG")) {
    EXPECT_EQ(
        m->stack[C][A][U][G] + m->stack[A][C][G][U] + m->stack[C][C][G][G] + m->hairpin["CCGAGG"],
        GetEnergy(kNNDBHairpin3));
  }

  EXPECT_EQ(m->stack[C][A][U][G] + m->stack[A][C][G][U] + m->stack[C][A][U][G] + m->au_penalty +
          m->terminal[A][C][C][U] + m->HairpinInitiation(6) + m->hairpin_all_c_a * 6 +
          m->hairpin_all_c_b,
      GetEnergy(kNNDBHairpin4));
  EXPECT_EQ(m->stack[C][G][C][G] + m->stack[G][G][C][C] + m->stack[G][G][U][C] + m->gu_penalty +
          m->terminal[G][G][G][U] + m->hairpin_gg_first_mismatch + m->HairpinInitiation(5) +
          m->hairpin_special_gu_closure,
      GetEnergy(kNNDBHairpin5));

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin1)), m);
    EXPECT_EQ(m->au_penalty + m->terminal[A][A][A][U] + m->HairpinInitiation(6), pc.Hairpin(3, 10));
  }

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin2)), m);
    EXPECT_EQ(m->au_penalty + m->terminal[A][G][G][U] + m->hairpin_gg_first_mismatch +
            m->HairpinInitiation(5),
        pc.Hairpin(3, 9));
  }

  if (m->hairpin.contains("CCGAGG")) {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin3)), m);
    EXPECT_EQ(m->hairpin["CCGAGG"], pc.Hairpin(3, 8));
  }

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin4)), m);
    EXPECT_EQ(m->au_penalty + m->terminal[A][C][C][U] + m->HairpinInitiation(6) +
            m->hairpin_all_c_a * 6 + m->hairpin_all_c_b,
        pc.Hairpin(3, 10));
  }

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin5)), m);
    EXPECT_EQ(m->gu_penalty + m->terminal[G][G][G][U] + m->hairpin_gg_first_mismatch +
            m->HairpinInitiation(5) + m->hairpin_special_gu_closure,
        pc.Hairpin(3, 9));
  }
}

TEST_P(EnergyTestBase, NNDBBulgeLoopExamples) {
  const auto& m = base_ms[GetParam()];

  EXPECT_EQ(m->stack[G][C][G][C] + m->stack[C][C][G][G] + m->BulgeInitiation(1) +
          m->bulge_special_c + m->stack[C][G][C][G] + m->HairpinInitiation(3) - E(R * T * log(3)),
      GetEnergy(kNNDBBulge1));

  EXPECT_EQ(m->stack[G][A][U][C] + m->au_penalty + m->BulgeInitiation(3) + m->HairpinInitiation(3),
      GetEnergy(kNNDBBulge2));
}

TEST_P(EnergyTestBase, NNDBInternalLoopExamples) {
  const auto& m = base_ms[GetParam()];

  EXPECT_EQ(m->stack[C][A][U][G] + m->stack[C][G][C][G] + m->InternalLoopInitiation(5) +
          std::min(m->internal_asym, NINIO_MAX_ASYM) + m->internal_2x3_mismatch[A][G][G][U] +
          m->internal_2x3_mismatch[G][G][A][C] + m->au_penalty + m->internal_au_penalty +
          m->HairpinInitiation(3),
      GetEnergy(kNNDBInternal2x3));
  EXPECT_EQ(m->stack[C][A][U][G] + m->stack[C][G][C][G] + m->internal_2x2[A][G][A][C][G][G][A][U] +
          m->au_penalty + m->HairpinInitiation(3),
      GetEnergy(kNNDBInternal2x2));
  EXPECT_EQ(m->stack[C][A][U][G] + m->stack[C][G][C][G] + m->InternalLoopInitiation(6) +
          std::min(4 * m->internal_asym, NINIO_MAX_ASYM) + m->au_penalty + m->internal_au_penalty +
          m->HairpinInitiation(3),
      GetEnergy(kNNDBInternal1x5));
}

TEST_P(EnergyTestBase, BaseCases) {
  const auto& m = base_ms[GetParam()];

  EXPECT_EQ(
      m->au_penalty + m->stack[G][A][U][C] + m->hairpin_init[3], GetEnergy("GAAAAUC", "((...))"));
  EXPECT_EQ(m->au_penalty + m->gu_penalty + m->stack[G][A][U][U] + m->hairpin_init[3],
      GetEnergy("GAAAAUU", "((...))"));
  EXPECT_EQ(m->au_penalty * 2 + m->HairpinInitiation(3) +
          std::min({ZERO_E, m->terminal[U][A][A][A],
              std::min(m->dangle3[U][A][A], m->dangle5[U][A][A])}),
      GetEnergy("AAAAAUA", ".(...)."));
  EXPECT_EQ(m->au_penalty * 2 + m->HairpinInitiation(3), GetEnergy("AAAAU", "(...)"));
  EXPECT_EQ(m->stack[G][C][G][C] + m->stack[C][U][A][G] + m->BulgeInitiation(1) +
          m->stack[U][G][C][A] + m->HairpinInitiation(3),
      GetEnergy(kBulge1));
  EXPECT_EQ(m->InternalLoopInitiation(5) + std::min(m->internal_asym, NINIO_MAX_ASYM) +
          m->internal_au_penalty + m->au_penalty * 2 + m->internal_2x3_mismatch[A][G][A][U] +
          m->internal_2x3_mismatch[C][A][A][G] + m->HairpinInitiation(3),
      GetEnergy(kInternal1));
}

INSTANTIATE_TEST_SUITE_P(EnergyModelTests, EnergyTestBase, testing::Range(0, NUM_TEST_MODELS));

#if ENERGY_PRECISION == 1

TEST(EnergyTestBase, T04) {
  auto m = base_t04;

  EXPECT_EQ(E(8.8), m->HairpinInitiation(87));
  EXPECT_EQ(E(6.8), m->BulgeInitiation(57));
  EXPECT_EQ(E(4.6), m->InternalLoopInitiation(67));

  const Precomp pc(Primary::FromSeq("GGGGAAACCCC"), m);
  EXPECT_EQ(E(-2.1 - 0.4 - 1.6), pc.min_mismatch_coax);
  EXPECT_EQ(E(-3.4), pc.min_flush_coax);
  EXPECT_EQ(E(-2.6), pc.min_twoloop_not_stack);

  Energy augubranch[4][4] = {{E(-0.6), E(-0.6), E(-0.6), E(0.5 - 0.6)},
      {E(-0.6), E(-0.6), E(-0.6), E(-0.6)}, {E(-0.6), E(-0.6), E(-0.6), E(0.5 - 0.6)},
      {E(0.5 - 0.6), E(-0.6), E(0.5 - 0.6), E(-0.6)}};
  EXPECT_EQ(sizeof(augubranch), sizeof(pc.augubranch));
  EXPECT_EQ(0, std::memcmp(augubranch, pc.augubranch, sizeof(augubranch)));
}

#elif ENERGY_PRECISION == 2

TEST(EnergyTestBase, T04) {
  auto m = base_t04;

  EXPECT_EQ(E(8.85), m->HairpinInitiation(87));
  EXPECT_EQ(E(6.79), m->BulgeInitiation(57));
  EXPECT_EQ(E(4.57), m->InternalLoopInitiation(67));

  const Precomp pc(Primary::FromSeq("GGGGAAACCCC"), m);
  EXPECT_EQ(E(-2.10 - 0.40 - 1.60), pc.min_mismatch_coax);
  EXPECT_EQ(E(-3.42), pc.min_flush_coax);
  EXPECT_EQ(E(-2.60), pc.min_twoloop_not_stack);

  Energy augubranch[4][4] = {{E(-0.60), E(-0.60), E(-0.60), E(0.50 - 0.60)},
      {E(-0.60), E(-0.60), E(-0.60), E(-0.60)}, {E(-0.60), E(-0.60), E(-0.60), E(0.50 - 0.60)},
      {E(0.50 - 0.60), E(-0.60), E(0.50 - 0.60), E(-0.60)}};
  EXPECT_EQ(sizeof(augubranch), sizeof(pc.augubranch));
  EXPECT_EQ(0, std::memcmp(augubranch, pc.augubranch, sizeof(augubranch)));
}

TEST(EnergyTestBase, T04Pseudofree) {
  auto m = base_t04;

  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "UGUCGAUACCCUGUCGAUA", "((((((((...))))))))");
    EXPECT_EQ(E(-3.88) + pseudofree, energy);
  }

  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "UAGGUCAGCCCCUGGUCUA", "((((((((...))))))))");
    EXPECT_EQ(E(-1.97) + pseudofree, energy);
  }

  // Penultimate stacking tests:

  // Hairpin + helices tests:
  // Test penultimate stacking not applied lonely pairs:
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGCCCCU", ".(...).");
    EXPECT_EQ(E(5.70) + pseudofree, energy);
  }

  // AU end on CG (0.44):
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGCCCCU", "((...))");
    EXPECT_EQ(E(5.32) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAGCCCCUU", ".((...)).");
    EXPECT_EQ(E(4.62) + pseudofree, energy);
  }

  // Counting twice, two AU on AUs (0.22):
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AACCCUU", "((...))");
    EXPECT_EQ(E(6.97) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAACCCUUU", ".((...)).");
    EXPECT_EQ(E(6.27) + pseudofree, energy);
  }

  // Special hairpin with lonely pair unaffected:
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "ACAGUGCU", "(......)");
    EXPECT_EQ(E(3.40) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAACAGUGCUUU", "..(......)..");
    EXPECT_EQ(E(2.70) + pseudofree, energy);
  }

  // Special hairpins:
  // GU end on AU at 0; AU end on GU at 1
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GACAGUGCUU", "((......))");
    EXPECT_EQ(E(2.13) + pseudofree, energy);
  }
  // GU end on AU at 1; AU end on GU at 2
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGACAGUGCUUU", ".((......)).");
    EXPECT_EQ(E(1.43) + pseudofree, energy);
  }

  // Single nuc bulge loop - treated as continuous.
  // AU end on GU at 0, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGCGAAAAUUUU", "((.((...))))");
    EXPECT_EQ(E(7.88) + pseudofree, energy);
  }
  // AU end on GU at 1, 5
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GAGCGAAAAUUUUU", ".((.((...)))).");
    EXPECT_EQ(E(7.18) + pseudofree, energy);
  }

  // Single nuc bulge, continuous single on outer:
  // GU end on GU at 0; AU end on GU at 3
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GCGAAAAUUU", "(.((...)))");
    EXPECT_EQ(E(8.43) + pseudofree, energy);
  }
  // GU end on GU at 1; AU end on GU at 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GGCGAAAAUUUU", ".(.((...))).");
    EXPECT_EQ(E(7.73) + pseudofree, energy);
  }

  // Single nuc bulge, continuous single on inner:
  // AU end on GU at 0, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGCAAAAUUU", "((.(...)))");
    EXPECT_EQ(E(8.38) + pseudofree, energy);
  }
  // AU end on GU at 1, 5
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GAGCAAAAUUUU", ".((.(...))).");
    EXPECT_EQ(E(7.68) + pseudofree, energy);
  }

  // Single nuc bulge, continuous single on both:
  // GU end on AU at 0; AU end on GU at 2
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GCAAAAUU", "(.(...))");
    EXPECT_EQ(E(8.93) + pseudofree, energy);
  }
  // GU end on AU at 1; AU end on GU at 3
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GGCAAAAUUU", ".(.(...)).");
    EXPECT_EQ(E(8.23) + pseudofree, energy);
  }

  // Multi nuc bulge loop:
  // AU end on GU at 0, 5; GU end on AU at 1, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGCCGAAAAUUUU", "((..((...))))");
    EXPECT_EQ(E(8.38) + pseudofree, energy);
  }
  // AU end on GU at 1, 6; GU end on AU at 2, 5
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GAGCCAGAAAUUUUU", ".((..((...)))).");
    EXPECT_EQ(E(8.40) + pseudofree, energy);
  }

  // Multi nuc bulge loop with lonely pairs:
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GCCGAAAUU", "(..(...))");
    EXPECT_EQ(E(10.20) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GGCCGAAAUUU", ".(..(...)).");
    EXPECT_EQ(E(9.50) + pseudofree, energy);
  }

  // Internal loop:
  // AU end on GU at 0, 4; GU end on AU at 1, 3
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGGGACCCUUCCCUU", "((.((...))...))");
    EXPECT_EQ(E(9.78) + pseudofree, energy);
  }
  // AU end on GU at 1, 5; GU end on AU at 2, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GAGGGACCCUUCCCUUU", ".((.((...))...)).");
    EXPECT_EQ(E(9.08) + pseudofree, energy);
  }

  // Internal loop with lonely pairs:
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAACCCUCCCU", "(.(...)...)");
    EXPECT_EQ(E(11.60) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAAACCCUCCCUU", ".(.(...)...).");
    EXPECT_EQ(E(10.90) + pseudofree, energy);
  }

  // Special 1x1 internal loop with helix affected by penultimate stack on external
  // sides only:
  // AU end on GU at 0, 4; GU end on AU at 1, 3
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGGGACCCUUCUU", "((.((...)).))");
    EXPECT_EQ(E(7.98) + pseudofree, energy);
  }
  // AU end on GU at 1, 5; GU end on AU at 2, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GAGGGACCCUUCUUU", ".((.((...)).)).");
    EXPECT_EQ(E(7.28) + pseudofree, energy);
  }

  // Special 1x1 internal loop with lonely pairs:
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAACCCUAU", "(.(...).)");
    EXPECT_EQ(E(9.80) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAAACCCUAUU", ".(.(...).).");
    EXPECT_EQ(E(9.10) + pseudofree, energy);
  }

  // Special 1x2 internal loop with helix affected by penultimate stack on external
  // sides only:
  // AU end on GU at 0, 4; GU end on AU at 1, 3
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGGGACCCUUCCUU", "((.((...))..))");
    EXPECT_EQ(E(9.78) + pseudofree, energy);
  }
  // AU end on GU at 1, 5; GU end on AU at 2, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GAGGGACCCUUCCUUU", ".((.((...))..)).");
    EXPECT_EQ(E(9.08) + pseudofree, energy);
  }

  // Special 1x2 internal loop with lonely pairs:
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAACCCUCAU", "(.(...)..)");
    EXPECT_EQ(E(11.60) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAAACCCUCAUU", ".(.(...)..).");
    EXPECT_EQ(E(10.90) + pseudofree, energy);
  }

  // Special 2x2 internal loop with helix affected by penultimate stack on external
  // sides only:
  // AU end on GU at 0, 5; GU end on AU at 1, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGGGGACCCUUCCUU", "((..((...))..))");
    EXPECT_EQ(E(9.18) + pseudofree, energy);
  }
  // AU end on GU at 1, 6; GU end on AU at 2, 5
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GAGGGGACCCUUCCUUU", ".((..((...))..)).");
    EXPECT_EQ(E(8.48) + pseudofree, energy);
  }

  // Special 2x2 internal loop with lonely pairs:
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAAACCCUAAU", "(..(...)..)");
    EXPECT_EQ(E(10.70) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAAAACCCUAAUU", ".(..(...)..).");
    EXPECT_EQ(E(10.00) + pseudofree, energy);
  }

  // Multiloop:
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "ACACCCUCCACCCUCCACCCUCU", "(.(...)..(...)..(...).)");
    EXPECT_EQ(E(27.70) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "AACACCCUCCACCCUCCACCCUCUU", ".(.(...)..(...)..(...).).");
    EXPECT_EQ(E(27.00) + pseudofree, energy);
  }

  // AU end on GU at 0; GU end on AU at 1
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "AGCACCCUCCACCCUCCACCCUCUU", "((.(...)..(...)..(...).))");
    EXPECT_EQ(E(27.15) + pseudofree, energy);
  }

  // AU end on GU at 1; GU end on AU at 2
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "AAGCACCCUCCACCCUCCACCCUCUUU", ".((.(...)..(...)..(...).)).");
    EXPECT_EQ(E(26.45) + pseudofree, energy);
  }

  // AU end on GU at 1, 4, 13; GU end on AU at 2, 5, 14
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "AAGCAGCCCUUCCAGCCCUUCCACCCUCUUU", ".((.((...))..((...))..(...).)).");
    EXPECT_EQ(E(25.35) + pseudofree, energy);
  }

  // Flush coax stack
  // GU end on AU at 0, 7; AU end on GU at 1, 8 (not taken as continuous)
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GACCCUUGACCCUU", "((...))((...))");
    EXPECT_EQ(E(13.26) + pseudofree, energy);
  }

  // GU end on GC at 1; GU end on AU at 7; AU end on GU at 8
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GGCCCUCGACCCUU", "((...))((...))");
    EXPECT_EQ(E(11.09) + pseudofree, energy);
  }

  // GU end on AU at 1, 8; AU end on GU at 2, 9
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGACCCUUGACCCUUU", ".((...))((...)).");
    EXPECT_EQ(E(12.86) + pseudofree, energy);
  }

  // GU end on GC at 2; GU end on AU at 8; AU end on GU at 9
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGGCCCUCGACCCUUU", ".((...))((...)).");
    EXPECT_EQ(E(11.09) + pseudofree, energy);
  }

  // Flush coax stack with lonely pairs:
  // Not counted as continuous, so no penultimate stacking.
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGCCCUACCCUU", ".(...)(...).");
    EXPECT_EQ(E(14.80) + pseudofree, energy);
  }

  // Mismatch mediated coax stack
  // GU end on AU at 1, 9; AU end on GU at 2, 10
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGACCCUUAGACCCUUU", ".((...)).((...)).");
    EXPECT_EQ(E(10.06) + pseudofree, energy);
  }

  // GU end on GC at 2; GU end on AU at 9; AU end on GU at 10
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGGCCCUCAGACCCUUU", ".((...)).((...)).");
    EXPECT_EQ(E(8.90) + pseudofree, energy);
  }

  // Mismatch mediated coax stack with lonely pairs:
  // Not counted as continuous, so no penultimate stacking.
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGCCCUAACCCUU", ".(...).(...).");
    EXPECT_EQ(E(12.60) + pseudofree, energy);
  }

  // Other tests:
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GCAAAGCC", "((...).)");
    EXPECT_EQ(E(4.45) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "CCCAAAAUG", ".(.(...))");
    EXPECT_EQ(E(5.71) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "UACAGA", "(....)");
    EXPECT_EQ(E(5.50) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGGGUCAUCCG", ".(((...))).");
    EXPECT_EQ(E(-0.59) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGAGAAACAAAU", "(..(...)...)");
    EXPECT_EQ(E(8.00) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "CGUUGCCUAAAAAGGAAACAAG", "(.............(...)..)");
    EXPECT_EQ(E(9.50) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "CCCGAAACAG", "(..(...).)");
    EXPECT_EQ(E(7.70) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GACAGAAACGCUGAAUC", "((..(...)......))");
    EXPECT_EQ(E(7.45) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "CUGAAACUGGAAACAGAAAUG", "(.(...)..(...).(...))");
    EXPECT_EQ(E(17.30) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))");
    EXPECT_EQ(E(18.25) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "AGCUAAAAACAAAGGUGAAACGU", "(..(...).(...)..(...).)");
    EXPECT_EQ(E(17.60) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "CUGAAACUGGAAACAGAAAUG", ".(.(...)(....)......)");
    EXPECT_EQ(E(13.10) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m,
        "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
        "(((((((((((.((...((((....))))..)).)))..((((..((((....))))...)))).))))))))....");
    EXPECT_EQ(E(-27.41) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "UCUGAGUAAAUUGCUACGCG", "(....)((...).......)");
    EXPECT_EQ(E(17.90) + pseudofree, energy);
  }

  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GGUCAAAGGUC", "((((...))))");
    EXPECT_EQ(E(3.63) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GGGGAAACCCC", "((((...))))");
    EXPECT_EQ(E(-4.38) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "UGACAAAGGCGA", "(..(...)...)");
    EXPECT_EQ(E(7.20) + pseudofree, energy);
  }
}

#endif

struct CtdTest {
  Primary r;
  Secondary s;
  Ctds ctd;
  BranchCtd branch_ctd;
  std::deque<int> branches;
};

namespace {

std::function<CtdTest(const Model::Ptr&)> CTD_TESTS[] = {
    [](const Model::Ptr&) -> CtdTest {
      return {.r = {}, .s = {}, .ctd = {}, .branch_ctd = {}, .branches = {}};
    },
    [](const Model::Ptr&) -> CtdTest {
      return {.r = Primary::FromSeq("A"),
          .s = Secondary::FromDb("."),
          .ctd = Ctds{CTD_NA},
          .branch_ctd = {},
          .branches = {}};
    },
    [](const Model::Ptr&) -> CtdTest {
      return {.r = Primary::FromSeq("AG"),
          .s = Secondary::FromDb(".."),
          .ctd = Ctds{CTD_NA, CTD_NA},
          .branch_ctd = {},
          .branches = {}};
    },
    [](const Model::Ptr&) -> CtdTest {
      return {.r = Primary::FromSeq("GUA"),
          .s = Secondary::FromDb("..."),
          .ctd = Ctds{CTD_NA, CTD_NA, CTD_NA},
          .branch_ctd = {},
          .branches = {}};
    },
    [](const Model::Ptr&) -> CtdTest {
      return {.r = Primary::FromSeq("GUAC"),
          .s = Secondary::FromDb("...."),
          .ctd = Ctds{CTD_NA, CTD_NA, CTD_NA, CTD_NA},
          .branch_ctd = {},
          .branches = {}};
    },
    // 3' dangle inside the branch.
    [](const Model::Ptr& m) -> CtdTest {
      return {.r = Primary::FromSeq("GAAAC"),
          .s = Secondary::FromDb("(...)"),
          .ctd = Ctds{CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_3_DANGLE},
          .branch_ctd = {{CTD_3_DANGLE, m->dangle3[G][A][C]}},
          .branches = {4}};
    },
    [](const Model::Ptr&) -> CtdTest {
      return {.r = Primary::FromSeq("GAAACAGAAAAUGGAAACCAGAAACA"),
          .s = Secondary::FromDb("(...).((...).(...)).(...)."),
          .ctd = Ctds(26),
          .branch_ctd = {},
          .branches = {}};
    },
    [](const Model::Ptr& m) -> CtdTest {
      return {.r = Primary::FromSeq("GAAACAGAAAAUGGAAACCAGAAACA"),
          .s = Secondary::FromDb("(...).((...).(...)).(...)."),
          .ctd = Ctds{CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_RC_WITH_NEXT, CTD_NA,
              CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
              CTD_NA, CTD_NA, CTD_RC_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA},
          .branch_ctd = {{CTD_UNUSED, ZERO_E}, {CTD_RC_WITH_NEXT, m->MismatchCoaxial(C, A, A, G)},
              {CTD_RC_WITH_PREV, m->MismatchCoaxial(C, A, A, G)}},
          .branches = {0, 6, 20}};
    },
    [](const Model::Ptr& m) -> CtdTest {
      return {.r = Primary::FromSeq("GAAACAGAAAAUGGAAACCAGAAACA"),
          .s = Secondary::FromDb("(...).((...).(...)).(...)."),
          .ctd = Ctds{CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_PREV,
              CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_5_DANGLE, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
              CTD_FCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA},
          .branch_ctd = {{CTD_FCOAX_WITH_NEXT, m->stack[G][A][U][C]},
              {CTD_FCOAX_WITH_PREV, m->stack[G][A][U][C]}, {CTD_5_DANGLE, m->dangle5[C][G][G]}},
          .branches = {18, 7, 13}};
    },
    [](const Model::Ptr& m) -> CtdTest {
      return {.r = Primary::FromSeq("GGAAACGAAACC"),
          .s = Secondary::FromDb("((...)(...))"),
          .ctd = Ctds{CTD_NA, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_NEXT,
              CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_PREV},
          .branch_ctd = {{CTD_UNUSED, ZERO_E}, {CTD_FCOAX_WITH_NEXT, m->stack[G][G][C][C]},
              {CTD_FCOAX_WITH_PREV, m->stack[G][G][C][C]}},
          .branches = {1, 6, 11}};
    },
    [](const Model::Ptr& m) -> CtdTest {
      return {.r = Primary::FromSeq("UUAGAAACGCAAAGAGGUCCAAAGA"),
          .s = Secondary::FromDb("(..(...).(...).....(...))"),
          .ctd = Ctds{CTD_NA, CTD_NA, CTD_NA, CTD_LCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
              CTD_NA, CTD_LCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
              CTD_NA, CTD_NA, CTD_FCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
              CTD_FCOAX_WITH_PREV},
          .branch_ctd = {{CTD_FCOAX_WITH_PREV, m->stack[U][C][G][A]},
              {CTD_LCOAX_WITH_NEXT, m->MismatchCoaxial(C, G, A, G)},
              {CTD_LCOAX_WITH_PREV, m->MismatchCoaxial(C, G, A, G)},
              {CTD_FCOAX_WITH_NEXT, m->stack[U][C][G][A]}},
          .branches = {24, 3, 9, 19}};
    }};
}  // namespace

class CtdsTestBase
    : public testing::TestWithParam<std::tuple<int, std::function<CtdTest(const Model::Ptr&)>>> {};

TEST_P(CtdsTestBase, BaseBranchBase) {
  const auto& m = base_ms[std::get<0>(GetParam())];
  auto ctd_test = std::get<1>(GetParam())(m);
  // Convert base representation to branch representation.
  BranchCtd computed_branch_ctd;
  auto computed_energy = AddBaseCtdsToBranchCtds(
      *m, ctd_test.r, ctd_test.s, ctd_test.ctd, ctd_test.branches, &computed_branch_ctd);
  Energy test_energy = ZERO_E;
  for (const auto& branch_ctd : ctd_test.branch_ctd) {
    // Make sure each branch energy is only represented once.
    if (branch_ctd.first == CTD_FCOAX_WITH_NEXT || branch_ctd.first == CTD_LCOAX_WITH_NEXT ||
        branch_ctd.first == CTD_RC_WITH_NEXT)
      continue;
    test_energy += branch_ctd.second;
  }
  EXPECT_EQ(test_energy, computed_energy);
  EXPECT_EQ(ctd_test.branch_ctd, computed_branch_ctd);
  // Convert back again and make sure it's the same.
  const Ctds prev_ctd = std::move(ctd_test.ctd);
  ctd_test.ctd.reset(prev_ctd.size());
  AddBranchCtdsToBaseCtds(ctd_test.branches, computed_branch_ctd, &ctd_test.ctd);
  EXPECT_EQ(prev_ctd, ctd_test.ctd);
}

INSTANTIATE_TEST_SUITE_P(CtdsTest, CtdsTestBase,
    testing::Combine(testing::Range(0, NUM_TEST_MODELS), testing::ValuesIn(CTD_TESTS)));

}  // namespace mrna::md::base
