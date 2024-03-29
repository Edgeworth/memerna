// Copyright 2016 E.
#include "model/energy.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>

#include "api/energy/energy.h"
#include "common_test.h"
#include "gtest/gtest.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "models/t22/energy/model.h"

namespace mrna::md::t22 {

Energy GetEnergy(const Model::Ptr& em, const std::tuple<Primary, Secondary>& s) {
  return em->TotalEnergy(std::get<Primary>(s), std::get<Secondary>(s), nullptr).energy;
}

Energy GetEnergy(const Model::Ptr& em, const std::string& r, const std::string& db) {
  return GetEnergy(em, {Primary::FromSeq(r), Secondary::FromDb(db)});
}

std::tuple<Energy, Energy> GetPseudofree(
    const Model::Ptr& base_em, const std::string& r, const std::string& db) {
  const auto paired_mul = E(100.0);
  const auto unpaired_mul = E(10.0);
  std::vector<Energy> pf_paired(r.size(), E(0.0));
  std::vector<Energy> pf_unpaired(r.size(), E(0.0));
  Energy extra_from_pseudofree = ZERO_E;
  for (int i = 0; i < int(r.size()); ++i) {
    pf_paired[i] = paired_mul * (i + 1);
    pf_unpaired[i] = unpaired_mul * (i + 1);
    if (db[i] == '.') {
      extra_from_pseudofree += pf_unpaired[i];
    } else {
      extra_from_pseudofree += pf_paired[i];
    }
  }

  auto em = base_em->CloneWithPseudofreeEnergy(pf_paired, pf_unpaired);
  auto energy = GetEnergy(em, {Primary::FromSeq(r), Secondary::FromDb(db)});
  return {energy, extra_from_pseudofree};
}

class T22ModelTest : public testing::TestWithParam<int> {
 public:
  static Energy GetEnergy(const std::string& r, const std::string& db) {
    return GetEnergy({Primary::FromSeq(r), Secondary::FromDb(db)});
  }

  static Energy GetEnergy(const std::tuple<Primary, Secondary>& s) {
    return test_t22_ems[GetParam()]
        ->TotalEnergy(std::get<Primary>(s), std::get<Secondary>(s), nullptr)
        .energy;
  }
};

TEST_P(T22ModelTest, MultiloopEnergy) {
  auto em = test_t22_ems[GetParam()];
  EXPECT_EQ(em->multiloop_hack_a + 4 * em->multiloop_hack_b, em->MultiloopInitiation(4));
}

TEST_P(T22ModelTest, NNDBHairpinLoopExamples) {
  auto em = test_t22_ems[GetParam()];

  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][A][U][G] + em->au_penalty +
          em->terminal[A][A][A][U] + em->HairpinInitiation(6) + em->penultimate_stack[C][A][U][G] +
          em->penultimate_stack[U][G][C][A],
      GetEnergy(kNNDBHairpin1));
  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][A][U][G] + em->au_penalty +
          em->terminal[A][G][G][U] + em->hairpin_gg_first_mismatch + em->HairpinInitiation(5) +
          em->penultimate_stack[C][A][U][G] + em->penultimate_stack[U][G][C][A],
      GetEnergy(kNNDBHairpin2));

  if (em->hairpin.contains("CCGAGG")) {
    EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][C][G][G] +
            em->hairpin["CCGAGG"] + em->penultimate_stack[C][C][G][G] +
            em->penultimate_stack[U][G][C][A],
        GetEnergy(kNNDBHairpin3));
  }

  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][A][U][G] + em->au_penalty +
          em->terminal[A][C][C][U] + em->HairpinInitiation(6) + em->hairpin_all_c_a * 6 +
          em->hairpin_all_c_b + em->penultimate_stack[C][A][U][G] +
          em->penultimate_stack[U][G][C][A],
      GetEnergy(kNNDBHairpin4));
  EXPECT_EQ(em->stack[C][G][C][G] + em->stack[G][G][C][C] + em->stack[G][G][U][C] + em->gu_penalty +
          em->terminal[G][G][G][U] + em->hairpin_gg_first_mismatch + em->HairpinInitiation(5) +
          em->hairpin_special_gu_closure + em->penultimate_stack[G][G][U][C] +
          em->penultimate_stack[C][G][C][G],
      GetEnergy(kNNDBHairpin5));
}

TEST_P(T22ModelTest, NNDBBulgeLoopExamples) {
  auto em = test_t22_ems[GetParam()];

  EXPECT_EQ(em->stack[G][C][G][C] + em->stack[C][C][G][G] + em->BulgeInitiation(1) +
          em->bulge_special_c + em->stack[C][G][C][G] + em->HairpinInitiation(3) -
          E(R * T * log(3)) + em->penultimate_stack[C][G][C][G] + em->penultimate_stack[G][C][G][C],
      GetEnergy(kNNDBBulge1));

  EXPECT_EQ(em->stack[G][A][U][C] + em->au_penalty + em->BulgeInitiation(3) +
          em->HairpinInitiation(3) + em->penultimate_stack[U][C][G][A] +
          em->penultimate_stack[G][A][U][C],
      GetEnergy(kNNDBBulge2));
}

TEST_P(T22ModelTest, NNDBInternalLoopExamples) {
  auto em = test_t22_ems[GetParam()];

  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[C][G][C][G] + em->InternalLoopInitiation(5) +
          std::min(em->internal_asym, NINIO_MAX_ASYM) + em->internal_2x3_mismatch[A][G][G][U] +
          em->internal_2x3_mismatch[G][G][A][C] + em->au_penalty + em->internal_au_penalty +
          em->HairpinInitiation(3) + em->penultimate_stack[C][G][C][G] +
          em->penultimate_stack[C][G][C][G] + em->penultimate_stack[C][A][U][G] +
          em->penultimate_stack[U][G][C][A],
      GetEnergy(kNNDBInternal2x3));
  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[C][G][C][G] +
          em->internal_2x2[A][G][A][C][G][G][A][U] + em->au_penalty + em->HairpinInitiation(3) +
          em->penultimate_stack[C][G][C][G] + em->penultimate_stack[C][G][C][G] +
          em->penultimate_stack[C][A][U][G] + em->penultimate_stack[U][G][C][A],
      GetEnergy(kNNDBInternal2x2));
  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[C][G][C][G] + em->InternalLoopInitiation(6) +
          std::min(4 * em->internal_asym, NINIO_MAX_ASYM) + em->au_penalty +
          em->internal_au_penalty + em->HairpinInitiation(3) + em->penultimate_stack[C][G][C][G] +
          em->penultimate_stack[C][G][C][G] + em->penultimate_stack[C][A][U][G] +
          em->penultimate_stack[U][G][C][A],
      GetEnergy(kNNDBInternal1x5));
}

TEST_P(T22ModelTest, BaseCases) {
  auto em = test_t22_ems[GetParam()];

  EXPECT_EQ(em->au_penalty + em->stack[G][A][U][C] + em->hairpin_init[3] +
          em->penultimate_stack[G][A][U][C] + em->penultimate_stack[U][C][G][A],
      GetEnergy("GAAAAUC", "((...))"));
  EXPECT_EQ(em->au_penalty + em->gu_penalty + em->stack[G][A][U][U] + em->hairpin_init[3] +
          em->penultimate_stack[G][A][U][U] + em->penultimate_stack[U][U][G][A],
      GetEnergy("GAAAAUU", "((...))"));
  EXPECT_EQ(em->au_penalty * 2 + em->HairpinInitiation(3) +
          std::min(ZERO_E,
              std::min(
                  em->terminal[U][A][A][A], std::min(em->dangle3[U][A][A], em->dangle5[U][A][A]))),
      GetEnergy("AAAAAUA", ".(...)."));
  EXPECT_EQ(em->au_penalty * 2 + em->HairpinInitiation(3), GetEnergy("AAAAU", "(...)"));
  EXPECT_EQ(em->stack[G][C][G][C] + em->stack[C][U][A][G] + em->BulgeInitiation(1) +
          em->stack[U][G][C][A] + em->HairpinInitiation(3) + em->penultimate_stack[U][G][C][A] +
          em->penultimate_stack[G][C][G][C],
      GetEnergy(kBulge1));
  EXPECT_EQ(em->InternalLoopInitiation(5) + std::min(em->internal_asym, NINIO_MAX_ASYM) +
          em->internal_au_penalty + em->au_penalty * 2 + em->internal_2x3_mismatch[A][G][A][U] +
          em->internal_2x3_mismatch[C][A][A][G] + em->HairpinInitiation(3),
      GetEnergy(kInternal1));
}

#if ENERGY_PRECISION == 2

TEST(T22P2ModelTest, T22P2) {
  auto em = t22p2;

  EXPECT_EQ(E(8.85), em->HairpinInitiation(87));
  EXPECT_EQ(E(6.79), em->BulgeInitiation(57));
  EXPECT_EQ(E(4.57), em->InternalLoopInitiation(67));

  // Example from https://doi.org/10.1093/nar/gkac261
  EXPECT_EQ(E(-3.65), GetEnergy(em, "UGUCGAUACCCUGUCGAUA", "((((((((...))))))))"));
  EXPECT_EQ(E(-4.05), GetEnergy(em, "UAGGUCAGCCCCUGGUCUA", "((((((((...))))))))"));

  // Penultimate stacking tests:

  // Hairpin + helices tests:
  // Test penultimate stacking not applied lonely pairs:
  EXPECT_EQ(E(5.70), GetEnergy(em, "AGCCCCU", ".(...)."));

  // AU end on CG (0.44):
  EXPECT_EQ(E(4.89 + 0.44), GetEnergy(em, "AGCCCCU", "((...))"));
  EXPECT_EQ(E(4.19 + 0.44), GetEnergy(em, "AAGCCCCUU", ".((...))."));

  // Counting twice, two AU on AUs (0.22):
  EXPECT_EQ(E(5.96 + 0.22 * 2), GetEnergy(em, "AACCCUU", "((...))"));
  EXPECT_EQ(E(5.26 + 0.22 * 2), GetEnergy(em, "AAACCCUUU", ".((...))."));

  // Special hairpin with lonely pair unaffected:
  EXPECT_EQ(E(2.40), GetEnergy(em, "ACAGUGCU", "(......)"));
  EXPECT_EQ(E(1.70), GetEnergy(em, "AAACAGUGCUUU", "..(......).."));

  // Special hairpins:
  // GU end on AU at 0; AU end on GU at 1
  EXPECT_EQ(E(1.80 - 0.31 - 0.71), GetEnergy(em, "GACAGUGCUU", "((......))"));
  // GU end on AU at 1; AU end on GU at 2
  EXPECT_EQ(E(1.10 - 0.31 - 0.71), GetEnergy(em, "AGACAGUGCUUU", ".((......))."));

  // Single nuc bulge loop - treated as continuous.
  // AU end on GU at 0, 4
  EXPECT_EQ(E(8.42 - 0.71 * 2), GetEnergy(em, "AGCGAAAAUUUU", "((.((...))))"));
  // AU end on GU at 1, 5
  EXPECT_EQ(E(7.72 - 0.71 * 2), GetEnergy(em, "GAGCGAAAAUUUUU", ".((.((...))))."));

  // Single nuc bulge, continuous single on outer:
  // GU end on GU at 0; AU end on GU at 3
  EXPECT_EQ(E(8.40 - 0.74 - 0.71), GetEnergy(em, "GCGAAAAUUU", "(.((...)))"));
  // GU end on GU at 1; AU end on GU at 4
  EXPECT_EQ(E(7.70 - 0.74 - 0.71), GetEnergy(em, "GGCGAAAAUUUU", ".(.((...)))."));

  // Single nuc bulge, continuous single on inner:
  // AU end on GU at 0, 4
  EXPECT_EQ(E(8.62 - 0.71 * 2), GetEnergy(em, "AGCAAAAUUU", "((.(...)))"));
  // AU end on GU at 1, 5
  EXPECT_EQ(E(7.92 - 0.71 * 2), GetEnergy(em, "GAGCAAAAUUUU", ".((.(...)))."));

  // Single nuc bulge, continuous single on both:
  // GU end on AU at 0; AU end on GU at 2
  EXPECT_EQ(E(8.60 - 0.31 - 0.71), GetEnergy(em, "GCAAAAUU", "(.(...))"));
  // GU end on AU at 1; AU end on GU at 3
  EXPECT_EQ(E(7.90 - 0.31 - 0.71), GetEnergy(em, "GGCAAAAUUU", ".(.(...))."));

  // Multi nuc bulge loop:
  // AU end on GU at 0, 5; GU end on AU at 1, 4
  EXPECT_EQ(E(7.62 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "AGCCGAAAAUUUU", "((..((...))))"));
  // AU end on GU at 1, 6; GU end on AU at 2, 5
  EXPECT_EQ(E(7.54 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "GAGCCAGAAAUUUUU", ".((..((...))))."));

  // Multi nuc bulge loop with lonely pairs:
  EXPECT_EQ(E(8.20), GetEnergy(em, "GCCGAAAUU", "(..(...))"));
  EXPECT_EQ(E(7.50), GetEnergy(em, "GGCCGAAAUUU", ".(..(...))."));

  // Internal loop:
  // AU end on GU at 0, 4; GU end on AU at 1, 3
  EXPECT_EQ(E(9.02 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "AGGGACCCUUCCCUU", "((.((...))...))"));
  // AU end on GU at 1, 5; GU end on AU at 2, 4
  EXPECT_EQ(E(8.32 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "GAGGGACCCUUCCCUUU", ".((.((...))...))."));

  // Internal loop with lonely pairs:
  EXPECT_EQ(E(9.60), GetEnergy(em, "AAACCCUCCCU", "(.(...)...)"));
  EXPECT_EQ(E(8.90), GetEnergy(em, "AAAACCCUCCCUU", ".(.(...)...)."));

  // Special 1x1 internal loop with helix affected by penultimate stack on external
  // sides only:
  // AU end on GU at 0, 4; GU end on AU at 1, 3
  EXPECT_EQ(E(7.22 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "AGGGACCCUUCUU", "((.((...)).))"));
  // AU end on GU at 1, 5; GU end on AU at 2, 4
  EXPECT_EQ(E(6.52 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "GAGGGACCCUUCUUU", ".((.((...)).))."));

  // Special 1x1 internal loop with lonely pairs:
  EXPECT_EQ(E(7.80), GetEnergy(em, "AAACCCUAU", "(.(...).)"));
  EXPECT_EQ(E(7.10), GetEnergy(em, "AAAACCCUAUU", ".(.(...).)."));

  // Special 1x2 internal loop with helix affected by penultimate stack on external
  // sides only:
  // AU end on GU at 0, 4; GU end on AU at 1, 3
  EXPECT_EQ(E(9.02 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "AGGGACCCUUCCUU", "((.((...))..))"));
  // AU end on GU at 1, 5; GU end on AU at 2, 4
  EXPECT_EQ(E(8.32 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "GAGGGACCCUUCCUUU", ".((.((...))..))."));

  // Special 1x2 internal loop with lonely pairs:
  EXPECT_EQ(E(9.60), GetEnergy(em, "AAACCCUCAU", "(.(...)..)"));
  EXPECT_EQ(E(8.90), GetEnergy(em, "AAAACCCUCAUU", ".(.(...)..)."));

  // Special 2x2 internal loop with helix affected by penultimate stack on external
  // sides only:
  // AU end on GU at 0, 5; GU end on AU at 1, 4
  EXPECT_EQ(E(8.42 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "AGGGGACCCUUCCUU", "((..((...))..))"));
  // AU end on GU at 1, 6; GU end on AU at 2, 5
  EXPECT_EQ(E(7.72 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "GAGGGGACCCUUCCUUU", ".((..((...))..))."));

  // Special 2x2 internal loop with lonely pairs:
  EXPECT_EQ(E(8.70), GetEnergy(em, "AAAACCCUAAU", "(..(...)..)"));
  EXPECT_EQ(E(8.00), GetEnergy(em, "AAAAACCCUAAUU", ".(..(...)..)."));

  // Multiloop:
  EXPECT_EQ(E(23.70), GetEnergy(em, "ACACCCUCCACCCUCCACCCUCU", "(.(...)..(...)..(...).)"));
  EXPECT_EQ(E(23.00), GetEnergy(em, "AACACCCUCCACCCUCCACCCUCUU", ".(.(...)..(...)..(...).)."));

  // AU end on GU at 0; GU end on AU at 1
  EXPECT_EQ(E(23.72 - 0.71 - 0.31),
      GetEnergy(em, "AGCACCCUCCACCCUCCACCCUCUU", "((.(...)..(...)..(...).))"));

  // AU end on GU at 1; GU end on AU at 2
  EXPECT_EQ(E(23.02 - 0.71 - 0.31),
      GetEnergy(em, "AAGCACCCUCCACCCUCCACCCUCUUU", ".((.(...)..(...)..(...).))."));

  // AU end on GU at 1, 4, 13; GU end on AU at 2, 5, 14
  EXPECT_EQ(E(23.06 - 0.71 * 3 - 0.31 * 3),
      GetEnergy(em, "AAGCAGCCCUUCCAGCCCUUCCACCCUCUUU", ".((.((...))..((...))..(...).))."));

  // Flush coax stack
  // GU end on AU at 0, 7; AU end on GU at 1, 8 (not taken as continuous)
  EXPECT_EQ(E(12.22 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "GACCCUUGACCCUU", "((...))((...))"));

  // GU end on GC at 1; GU end on AU at 7; AU end on GU at 8
  EXPECT_EQ(E(10.35 + 0.13 - 0.71 - 0.31), GetEnergy(em, "GGCCCUCGACCCUU", "((...))((...))"));

  // GU end on AU at 1, 8; AU end on GU at 2, 9
  EXPECT_EQ(E(12.20 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "AGACCCUUGACCCUUU", ".((...))((...))."));

  // GU end on GC at 2; GU end on AU at 8; AU end on GU at 9
  EXPECT_EQ(E(10.35 + 0.13 - 0.71 - 0.31), GetEnergy(em, "AGGCCCUCGACCCUUU", ".((...))((...))."));

  // Flush coax stack with lonely pairs:
  // Not counted as continuous, so no penultimate stacking.
  EXPECT_EQ(E(13.40), GetEnergy(em, "AGCCCUACCCUU", ".(...)(...)."));

  // Mismatch mediated coax stack
  // GU end on AU at 1, 9; AU end on GU at 2, 10
  EXPECT_EQ(E(9.40 - 0.71 * 2 - 0.31 * 2), GetEnergy(em, "AGACCCUUAGACCCUUU", ".((...)).((...))."));

  // GU end on GC at 2; GU end on AU at 9; AU end on GU at 10
  EXPECT_EQ(E(7.80 + 0.13 - 0.71 - 0.31), GetEnergy(em, "AGGCCCUCAGACCCUUU", ".((...)).((...))."));

  // Mismatch mediated coax stack with lonely pairs:
  // Not counted as continuous, so no penultimate stacking.
  EXPECT_EQ(E(10.60), GetEnergy(em, "AGCCCUAACCCUU", ".(...).(...)."));

  // Other tests:
  EXPECT_EQ(E(4.41), GetEnergy(em, "GCAAAGCC", "((...).)"));
  EXPECT_EQ(E(5.69), GetEnergy(em, "CCCAAAAUG", ".(.(...))"));
  EXPECT_EQ(E(4.50), GetEnergy(em, "UACAGA", "(....)"));
  EXPECT_EQ(E(-1.25), GetEnergy(em, "AGGGUCAUCCG", ".(((...)))."));
  EXPECT_EQ(E(7.00), GetEnergy(em, "AGAGAAACAAAU", "(..(...)...)"));
  EXPECT_EQ(E(9.50), GetEnergy(em, "CGUUGCCUAAAAAGGAAACAAG", "(.............(...)..)"));
  EXPECT_EQ(E(7.70), GetEnergy(em, "CCCGAAACAG", "(..(...).)"));
  EXPECT_EQ(E(7.32), GetEnergy(em, "GACAGAAACGCUGAAUC", "((..(...)......))"));
  EXPECT_EQ(E(16.30), GetEnergy(em, "CUGAAACUGGAAACAGAAAUG", "(.(...)..(...).(...))"));
  EXPECT_EQ(E(17.18), GetEnergy(em, "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))"));
  EXPECT_EQ(E(15.60), GetEnergy(em, "AGCUAAAAACAAAGGUGAAACGU", "(..(...).(...)..(...).)"));
  EXPECT_EQ(E(11.10), GetEnergy(em, "CUGAAACUGGAAACAGAAAUG", ".(.(...)(....)......)"));
  EXPECT_EQ(E(-29.04),
      GetEnergy(em, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
          "(((((((((((.((...((((....))))..)).)))..((((..((((....))))...)))).))))))))...."));
  EXPECT_EQ(E(14.71), GetEnergy(em, "UCUGAGUAAAUUGCUACGCG", "(....)((...).......)"));
  EXPECT_EQ(E(-45.38), GetEnergy(em, k16sHSapiens3));

  EXPECT_EQ(E(1.61), GetEnergy(em, "GGUCAAAGGUC", "((((...))))"));
  EXPECT_EQ(E(-4.44), GetEnergy(em, "GGGGAAACCCC", "((((...))))"));
  EXPECT_EQ(E(6.20), GetEnergy(em, "UGACAAAGGCGA", "(..(...)...)"));

  // NNDB flush coax
  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][A][U][G] +
          2 * em->au_penalty + 2 * em->HairpinInitiation(3) + em->penultimate_stack[U][G][C][A] +
          em->penultimate_stack[C][A][U][G] + em->penultimate_stack[C][A][U][G] +
          em->penultimate_stack[U][G][C][A],
      GetEnergy(em, "GUGAAACACAAAAUGA", ".((...))((...))."));

  // NNDB T99 Multiloop example
  EXPECT_EQ(em->stack[G][A][U][C] + em->terminal[C][G][A][G] + em->coax_mismatch_non_contiguous +
          3 * em->HairpinInitiation(3) + em->MultiloopInitiation(4) + 2 * em->au_penalty,
      GetEnergy(em, "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))"));
}

TEST(T22P2ModelTest, T22P2Pseudofree) {
  auto em = t22p2;

  {
    // Example from https://doi.org/10.1093/nar/gkac261
    const auto& [energy, pseudofree] =
        GetPseudofree(em, "UGUCGAUACCCUGUCGAUA", "((((((((...))))))))");
    EXPECT_EQ(E(-3.65) + pseudofree, energy);
  }

  {
    const auto& [energy, pseudofree] =
        GetPseudofree(em, "UAGGUCAGCCCCUGGUCUA", "((((((((...))))))))");
    EXPECT_EQ(E(-4.05) + pseudofree, energy);
  }

  // Penultimate stacking tests:

  // Hairpin + helices tests:
  // Test penultimate stacking not applied lonely pairs:
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AGCCCCU", ".(...).");
    EXPECT_EQ(E(5.70) + pseudofree, energy);
  }

  // AU end on CG (0.44):
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AGCCCCU", "((...))");
    EXPECT_EQ(E(4.89 + 0.44) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AAGCCCCUU", ".((...)).");
    EXPECT_EQ(E(4.19 + 0.44) + pseudofree, energy);
  }

  // Counting twice, two AU on AUs (0.22):
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AACCCUU", "((...))");
    EXPECT_EQ(E(5.96 + 0.22 * 2) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AAACCCUUU", ".((...)).");
    EXPECT_EQ(E(5.26 + 0.22 * 2) + pseudofree, energy);
  }

  // Special hairpin with lonely pair unaffected:
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "ACAGUGCU", "(......)");
    EXPECT_EQ(E(2.40) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AAACAGUGCUUU", "..(......)..");
    EXPECT_EQ(E(1.70) + pseudofree, energy);
  }

  // Special hairpins:
  // GU end on AU at 0; AU end on GU at 1
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GACAGUGCUU", "((......))");
    EXPECT_EQ(E(1.80 - 0.31 - 0.71) + pseudofree, energy);
  }
  // GU end on AU at 1; AU end on GU at 2
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AGACAGUGCUUU", ".((......)).");
    EXPECT_EQ(E(1.10 - 0.31 - 0.71) + pseudofree, energy);
  }

  // Single nuc bulge loop - treated as continuous.
  // AU end on GU at 0, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AGCGAAAAUUUU", "((.((...))))");
    EXPECT_EQ(E(8.42 - 0.71 * 2) + pseudofree, energy);
  }
  // AU end on GU at 1, 5
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GAGCGAAAAUUUUU", ".((.((...)))).");
    EXPECT_EQ(E(7.72 - 0.71 * 2) + pseudofree, energy);
  }

  // Single nuc bulge, continuous single on outer:
  // GU end on GU at 0; AU end on GU at 3
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GCGAAAAUUU", "(.((...)))");
    EXPECT_EQ(E(8.40 - 0.74 - 0.71) + pseudofree, energy);
  }
  // GU end on GU at 1; AU end on GU at 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GGCGAAAAUUUU", ".(.((...))).");
    EXPECT_EQ(E(7.70 - 0.74 - 0.71) + pseudofree, energy);
  }

  // Single nuc bulge, continuous single on inner:
  // AU end on GU at 0, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AGCAAAAUUU", "((.(...)))");
    EXPECT_EQ(E(8.62 - 0.71 * 2) + pseudofree, energy);
  }
  // AU end on GU at 1, 5
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GAGCAAAAUUUU", ".((.(...))).");
    EXPECT_EQ(E(7.92 - 0.71 * 2) + pseudofree, energy);
  }

  // Single nuc bulge, continuous single on both:
  // GU end on AU at 0; AU end on GU at 2
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GCAAAAUU", "(.(...))");
    EXPECT_EQ(E(8.60 - 0.31 - 0.71) + pseudofree, energy);
  }
  // GU end on AU at 1; AU end on GU at 3
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GGCAAAAUUU", ".(.(...)).");
    EXPECT_EQ(E(7.90 - 0.31 - 0.71) + pseudofree, energy);
  }

  // Multi nuc bulge loop:
  // AU end on GU at 0, 5; GU end on AU at 1, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AGCCGAAAAUUUU", "((..((...))))");
    EXPECT_EQ(E(7.62 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }
  // AU end on GU at 1, 6; GU end on AU at 2, 5
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GAGCCAGAAAUUUUU", ".((..((...)))).");
    EXPECT_EQ(E(7.54 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }

  // Multi nuc bulge loop with lonely pairs:
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GCCGAAAUU", "(..(...))");
    EXPECT_EQ(E(8.20) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GGCCGAAAUUU", ".(..(...)).");
    EXPECT_EQ(E(7.50) + pseudofree, energy);
  }

  // Internal loop:
  // AU end on GU at 0, 4; GU end on AU at 1, 3
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AGGGACCCUUCCCUU", "((.((...))...))");
    EXPECT_EQ(E(9.02 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }
  // AU end on GU at 1, 5; GU end on AU at 2, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GAGGGACCCUUCCCUUU", ".((.((...))...)).");
    EXPECT_EQ(E(8.32 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }

  // Internal loop with lonely pairs:
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AAACCCUCCCU", "(.(...)...)");
    EXPECT_EQ(E(9.60) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AAAACCCUCCCUU", ".(.(...)...).");
    EXPECT_EQ(E(8.90) + pseudofree, energy);
  }

  // Special 1x1 internal loop with helix affected by penultimate stack on external
  // sides only:
  // AU end on GU at 0, 4; GU end on AU at 1, 3
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AGGGACCCUUCUU", "((.((...)).))");
    EXPECT_EQ(E(7.22 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }
  // AU end on GU at 1, 5; GU end on AU at 2, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GAGGGACCCUUCUUU", ".((.((...)).)).");
    EXPECT_EQ(E(6.52 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }

  // Special 1x1 internal loop with lonely pairs:
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AAACCCUAU", "(.(...).)");
    EXPECT_EQ(E(7.80) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AAAACCCUAUU", ".(.(...).).");
    EXPECT_EQ(E(7.10) + pseudofree, energy);
  }

  // Special 1x2 internal loop with helix affected by penultimate stack on external
  // sides only:
  // AU end on GU at 0, 4; GU end on AU at 1, 3
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AGGGACCCUUCCUU", "((.((...))..))");
    EXPECT_EQ(E(9.02 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }
  // AU end on GU at 1, 5; GU end on AU at 2, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GAGGGACCCUUCCUUU", ".((.((...))..)).");
    EXPECT_EQ(E(8.32 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }

  // Special 1x2 internal loop with lonely pairs:
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AAACCCUCAU", "(.(...)..)");
    EXPECT_EQ(E(9.60) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AAAACCCUCAUU", ".(.(...)..).");
    EXPECT_EQ(E(8.90) + pseudofree, energy);
  }

  // Special 2x2 internal loop with helix affected by penultimate stack on external
  // sides only:
  // AU end on GU at 0, 5; GU end on AU at 1, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AGGGGACCCUUCCUU", "((..((...))..))");
    EXPECT_EQ(E(8.42 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }
  // AU end on GU at 1, 6; GU end on AU at 2, 5
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GAGGGGACCCUUCCUUU", ".((..((...))..)).");
    EXPECT_EQ(E(7.72 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }

  // Special 2x2 internal loop with lonely pairs:
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AAAACCCUAAU", "(..(...)..)");
    EXPECT_EQ(E(8.70) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AAAAACCCUAAUU", ".(..(...)..).");
    EXPECT_EQ(E(8.00) + pseudofree, energy);
  }

  // Multiloop:
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(em, "ACACCCUCCACCCUCCACCCUCU", "(.(...)..(...)..(...).)");
    EXPECT_EQ(E(23.70) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(em, "AACACCCUCCACCCUCCACCCUCUU", ".(.(...)..(...)..(...).).");
    EXPECT_EQ(E(23.00) + pseudofree, energy);
  }

  // AU end on GU at 0; GU end on AU at 1
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(em, "AGCACCCUCCACCCUCCACCCUCUU", "((.(...)..(...)..(...).))");
    EXPECT_EQ(E(23.72 - 0.71 - 0.31) + pseudofree, energy);
  }

  // AU end on GU at 1; GU end on AU at 2
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(em, "AAGCACCCUCCACCCUCCACCCUCUUU", ".((.(...)..(...)..(...).)).");
    EXPECT_EQ(E(23.02 - 0.71 - 0.31) + pseudofree, energy);
  }

  // AU end on GU at 1, 4, 13; GU end on AU at 2, 5, 14
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(em, "AAGCAGCCCUUCCAGCCCUUCCACCCUCUUU", ".((.((...))..((...))..(...).)).");
    EXPECT_EQ(E(23.06 - 0.71 * 3 - 0.31 * 3) + pseudofree, energy);
  }

  // Flush coax stack
  // GU end on AU at 0, 7; AU end on GU at 1, 8 (not taken as continuous)
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GACCCUUGACCCUU", "((...))((...))");
    EXPECT_EQ(E(12.22 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }

  // GU end on GC at 1; GU end on AU at 7; AU end on GU at 8
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GGCCCUCGACCCUU", "((...))((...))");
    EXPECT_EQ(E(10.35 + 0.13 - 0.71 - 0.31) + pseudofree, energy);
  }

  // GU end on AU at 1, 8; AU end on GU at 2, 9
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AGACCCUUGACCCUUU", ".((...))((...)).");
    EXPECT_EQ(E(12.20 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }

  // GU end on GC at 2; GU end on AU at 8; AU end on GU at 9
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AGGCCCUCGACCCUUU", ".((...))((...)).");
    EXPECT_EQ(E(10.35 + 0.13 - 0.71 - 0.31) + pseudofree, energy);
  }

  // Flush coax stack with lonely pairs:
  // Not counted as continuous, so no penultimate stacking.
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AGCCCUACCCUU", ".(...)(...).");
    EXPECT_EQ(E(13.40) + pseudofree, energy);
  }

  // Mismatch mediated coax stack
  // GU end on AU at 1, 9; AU end on GU at 2, 10
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AGACCCUUAGACCCUUU", ".((...)).((...)).");
    EXPECT_EQ(E(9.40 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }

  // GU end on GC at 2; GU end on AU at 9; AU end on GU at 10
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AGGCCCUCAGACCCUUU", ".((...)).((...)).");
    EXPECT_EQ(E(7.80 + 0.13 - 0.71 - 0.31) + pseudofree, energy);
  }

  // Mismatch mediated coax stack with lonely pairs:
  // Not counted as continuous, so no penultimate stacking.
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AGCCCUAACCCUU", ".(...).(...).");
    EXPECT_EQ(E(10.60) + pseudofree, energy);
  }

  // Other tests:
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GCAAAGCC", "((...).)");
    EXPECT_EQ(E(4.41) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "CCCAAAAUG", ".(.(...))");
    EXPECT_EQ(E(5.69) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "UACAGA", "(....)");
    EXPECT_EQ(E(4.50) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AGGGUCAUCCG", ".(((...))).");
    EXPECT_EQ(E(-1.25) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "AGAGAAACAAAU", "(..(...)...)");
    EXPECT_EQ(E(7.00) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(em, "CGUUGCCUAAAAAGGAAACAAG", "(.............(...)..)");
    EXPECT_EQ(E(9.50) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "CCCGAAACAG", "(..(...).)");
    EXPECT_EQ(E(7.70) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GACAGAAACGCUGAAUC", "((..(...)......))");
    EXPECT_EQ(E(7.32) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(em, "CUGAAACUGGAAACAGAAAUG", "(.(...)..(...).(...))");
    EXPECT_EQ(E(16.30) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(em, "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))");
    EXPECT_EQ(E(17.18) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(em, "AGCUAAAAACAAAGGUGAAACGU", "(..(...).(...)..(...).)");
    EXPECT_EQ(E(15.60) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(em, "CUGAAACUGGAAACAGAAAUG", ".(.(...)(....)......)");
    EXPECT_EQ(E(11.10) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(em,
        "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
        "(((((((((((.((...((((....))))..)).)))..((((..((((....))))...)))).))))))))....");
    EXPECT_EQ(E(-29.04) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(em, "UCUGAGUAAAUUGCUACGCG", "(....)((...).......)");
    EXPECT_EQ(E(14.71) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(
        em, std::get<Primary>(k16sHSapiens3).ToSeq(), std::get<Secondary>(k16sHSapiens3).ToDb());
    EXPECT_EQ(E(-45.38) + pseudofree, energy);
  }

  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GGUCAAAGGUC", "((((...))))");
    EXPECT_EQ(E(1.61) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GGGGAAACCCC", "((((...))))");
    EXPECT_EQ(E(-4.44) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "UGACAAAGGCGA", "(..(...)...)");
    EXPECT_EQ(E(6.20) + pseudofree, energy);
  }

  // NNDB flush coax
  {
    const auto& [energy, pseudofree] = GetPseudofree(em, "GUGAAACACAAAAUGA", ".((...))((...)).");
    EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][A][U][G] +
            2 * em->au_penalty + 2 * em->HairpinInitiation(3) + em->penultimate_stack[U][G][C][A] +
            em->penultimate_stack[C][A][U][G] + em->penultimate_stack[C][A][U][G] +
            em->penultimate_stack[U][G][C][A] + pseudofree,
        energy);
  }

  // NNDB T99 Multiloop example
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(em, "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))");
    EXPECT_EQ(em->stack[G][A][U][C] + em->terminal[C][G][A][G] + em->coax_mismatch_non_contiguous +
            3 * em->HairpinInitiation(3) + em->MultiloopInitiation(4) + 2 * em->au_penalty +
            pseudofree,
        energy);
  }
}

#endif

INSTANTIATE_TEST_SUITE_P(EnergyModelTests, T22ModelTest, testing::Range(0, NUM_TEST_MODELS));

}  // namespace mrna::md::t22
