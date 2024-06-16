// Copyright 2024 Eliot Courtney.
#include "model/energy.h"

#include <cmath>
#include <string>

#include "gtest/gtest.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "tests/init.h"
#include "tests/util.h"

namespace mrna {

class EnergyTestT22 : public testing::TestWithParam<int> {};

#if ENERGY_PRECISION == 1

GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(EnergyTestT22);

#elif ENERGY_PRECISION == 2

TEST_P(EnergyTestT22, T22P2) {
  auto m = t22_ms[GetParam()];

  // Example from https://doi.org/10.1093/nar/gkac261
  EXPECT_EQ(E(-3.65), GetEnergy(m, "UGUCGAUACCCUGUCGAUA", "((((((((...))))))))"));
  EXPECT_EQ(E(-4.05), GetEnergy(m, "UAGGUCAGCCCCUGGUCUA", "((((((((...))))))))"));

  // Penultimate stacking tests:

  // Hairpin + helices tests:
  // Test penultimate stacking not applied lonely pairs:
  EXPECT_EQ(E(5.70), GetEnergy(m, "AGCCCCU", ".(...)."));

  // AU end on CG (0.44):
  EXPECT_EQ(E(4.89 + 0.44), GetEnergy(m, "AGCCCCU", "((...))"));
  EXPECT_EQ(E(4.19 + 0.44), GetEnergy(m, "AAGCCCCUU", ".((...))."));

  // Counting twice, two AU on AUs (0.22):
  EXPECT_EQ(E(5.96 + 0.22 * 2), GetEnergy(m, "AACCCUU", "((...))"));
  EXPECT_EQ(E(5.26 + 0.22 * 2), GetEnergy(m, "AAACCCUUU", ".((...))."));

  // Special hairpin with lonely pair unaffected:
  EXPECT_EQ(E(2.40), GetEnergy(m, "ACAGUGCU", "(......)"));
  EXPECT_EQ(E(1.70), GetEnergy(m, "AAACAGUGCUUU", "..(......).."));

  // Special hairpins:
  // GU end on AU at 0; AU end on GU at 1
  EXPECT_EQ(E(1.80 - 0.31 - 0.71), GetEnergy(m, "GACAGUGCUU", "((......))"));
  // GU end on AU at 1; AU end on GU at 2
  EXPECT_EQ(E(1.10 - 0.31 - 0.71), GetEnergy(m, "AGACAGUGCUUU", ".((......))."));

  // Single nuc bulge loop - treated as continuous.
  // AU end on GU at 0, 4
  EXPECT_EQ(E(8.42 - 0.71 * 2), GetEnergy(m, "AGCGAAAAUUUU", "((.((...))))"));
  // AU end on GU at 1, 5
  EXPECT_EQ(E(7.72 - 0.71 * 2), GetEnergy(m, "GAGCGAAAAUUUUU", ".((.((...))))."));

  // Single nuc bulge, continuous single on outer:
  // GU end on GU at 0; AU end on GU at 3
  EXPECT_EQ(E(8.40 - 0.74 - 0.71), GetEnergy(m, "GCGAAAAUUU", "(.((...)))"));
  // GU end on GU at 1; AU end on GU at 4
  EXPECT_EQ(E(7.70 - 0.74 - 0.71), GetEnergy(m, "GGCGAAAAUUUU", ".(.((...)))."));

  // Single nuc bulge, continuous single on inner:
  // AU end on GU at 0, 4
  EXPECT_EQ(E(8.62 - 0.71 * 2), GetEnergy(m, "AGCAAAAUUU", "((.(...)))"));
  // AU end on GU at 1, 5
  EXPECT_EQ(E(7.92 - 0.71 * 2), GetEnergy(m, "GAGCAAAAUUUU", ".((.(...)))."));

  // Single nuc bulge, continuous single on both:
  // GU end on AU at 0; AU end on GU at 2
  EXPECT_EQ(E(8.60 - 0.31 - 0.71), GetEnergy(m, "GCAAAAUU", "(.(...))"));
  // GU end on AU at 1; AU end on GU at 3
  EXPECT_EQ(E(7.90 - 0.31 - 0.71), GetEnergy(m, "GGCAAAAUUU", ".(.(...))."));

  // Multi nuc bulge loop:
  // AU end on GU at 0, 5; GU end on AU at 1, 4
  EXPECT_EQ(E(7.62 - 0.71 * 2 - 0.31 * 2), GetEnergy(m, "AGCCGAAAAUUUU", "((..((...))))"));
  // AU end on GU at 1, 6; GU end on AU at 2, 5
  EXPECT_EQ(E(7.54 - 0.71 * 2 - 0.31 * 2), GetEnergy(m, "GAGCCAGAAAUUUUU", ".((..((...))))."));

  // Multi nuc bulge loop with lonely pairs:
  EXPECT_EQ(E(8.20), GetEnergy(m, "GCCGAAAUU", "(..(...))"));
  EXPECT_EQ(E(7.50), GetEnergy(m, "GGCCGAAAUUU", ".(..(...))."));

  // Internal loop:
  // AU end on GU at 0, 4; GU end on AU at 1, 3
  EXPECT_EQ(E(9.02 - 0.71 * 2 - 0.31 * 2), GetEnergy(m, "AGGGACCCUUCCCUU", "((.((...))...))"));
  // AU end on GU at 1, 5; GU end on AU at 2, 4
  EXPECT_EQ(E(8.32 - 0.71 * 2 - 0.31 * 2), GetEnergy(m, "GAGGGACCCUUCCCUUU", ".((.((...))...))."));

  // Internal loop with lonely pairs:
  EXPECT_EQ(E(9.60), GetEnergy(m, "AAACCCUCCCU", "(.(...)...)"));
  EXPECT_EQ(E(8.90), GetEnergy(m, "AAAACCCUCCCUU", ".(.(...)...)."));

  // Special 1x1 internal loop with helix affected by penultimate stack on external
  // sides only:
  // AU end on GU at 0, 4; GU end on AU at 1, 3
  EXPECT_EQ(E(7.22 - 0.71 * 2 - 0.31 * 2), GetEnergy(m, "AGGGACCCUUCUU", "((.((...)).))"));
  // AU end on GU at 1, 5; GU end on AU at 2, 4
  EXPECT_EQ(E(6.52 - 0.71 * 2 - 0.31 * 2), GetEnergy(m, "GAGGGACCCUUCUUU", ".((.((...)).))."));

  // Special 1x1 internal loop with lonely pairs:
  EXPECT_EQ(E(7.80), GetEnergy(m, "AAACCCUAU", "(.(...).)"));
  EXPECT_EQ(E(7.10), GetEnergy(m, "AAAACCCUAUU", ".(.(...).)."));

  // Special 1x2 internal loop with helix affected by penultimate stack on external
  // sides only:
  // AU end on GU at 0, 4; GU end on AU at 1, 3
  EXPECT_EQ(E(9.02 - 0.71 * 2 - 0.31 * 2), GetEnergy(m, "AGGGACCCUUCCUU", "((.((...))..))"));
  // AU end on GU at 1, 5; GU end on AU at 2, 4
  EXPECT_EQ(E(8.32 - 0.71 * 2 - 0.31 * 2), GetEnergy(m, "GAGGGACCCUUCCUUU", ".((.((...))..))."));

  // Special 1x2 internal loop with lonely pairs:
  EXPECT_EQ(E(9.60), GetEnergy(m, "AAACCCUCAU", "(.(...)..)"));
  EXPECT_EQ(E(8.90), GetEnergy(m, "AAAACCCUCAUU", ".(.(...)..)."));

  // Special 2x2 internal loop with helix affected by penultimate stack on external
  // sides only:
  // AU end on GU at 0, 5; GU end on AU at 1, 4
  EXPECT_EQ(E(8.42 - 0.71 * 2 - 0.31 * 2), GetEnergy(m, "AGGGGACCCUUCCUU", "((..((...))..))"));
  // AU end on GU at 1, 6; GU end on AU at 2, 5
  EXPECT_EQ(E(7.72 - 0.71 * 2 - 0.31 * 2), GetEnergy(m, "GAGGGGACCCUUCCUUU", ".((..((...))..))."));

  // Special 2x2 internal loop with lonely pairs:
  EXPECT_EQ(E(8.70), GetEnergy(m, "AAAACCCUAAU", "(..(...)..)"));
  EXPECT_EQ(E(8.00), GetEnergy(m, "AAAAACCCUAAUU", ".(..(...)..)."));

  // Multiloop:
  EXPECT_EQ(E(23.70), GetEnergy(m, "ACACCCUCCACCCUCCACCCUCU", "(.(...)..(...)..(...).)"));
  EXPECT_EQ(E(23.00), GetEnergy(m, "AACACCCUCCACCCUCCACCCUCUU", ".(.(...)..(...)..(...).)."));

  // AU end on GU at 0; GU end on AU at 1
  EXPECT_EQ(E(23.72 - 0.71 - 0.31),
      GetEnergy(m, "AGCACCCUCCACCCUCCACCCUCUU", "((.(...)..(...)..(...).))"));

  // AU end on GU at 1; GU end on AU at 2
  EXPECT_EQ(E(23.02 - 0.71 - 0.31),
      GetEnergy(m, "AAGCACCCUCCACCCUCCACCCUCUUU", ".((.(...)..(...)..(...).))."));

  // AU end on GU at 1, 4, 13; GU end on AU at 2, 5, 14
  EXPECT_EQ(E(23.06 - 0.71 * 3 - 0.31 * 3),
      GetEnergy(m, "AAGCAGCCCUUCCAGCCCUUCCACCCUCUUU", ".((.((...))..((...))..(...).))."));

  // Flush coax stack
  // GU end on AU at 0, 7; AU end on GU at 1, 8 (not taken as continuous)
  EXPECT_EQ(E(12.22 - 0.71 * 2 - 0.31 * 2), GetEnergy(m, "GACCCUUGACCCUU", "((...))((...))"));

  // GU end on GC at 1; GU end on AU at 7; AU end on GU at 8
  EXPECT_EQ(E(10.35 + 0.13 - 0.71 - 0.31), GetEnergy(m, "GGCCCUCGACCCUU", "((...))((...))"));

  // GU end on AU at 1, 8; AU end on GU at 2, 9
  EXPECT_EQ(E(12.20 - 0.71 * 2 - 0.31 * 2), GetEnergy(m, "AGACCCUUGACCCUUU", ".((...))((...))."));

  // GU end on GC at 2; GU end on AU at 8; AU end on GU at 9
  EXPECT_EQ(E(10.35 + 0.13 - 0.71 - 0.31), GetEnergy(m, "AGGCCCUCGACCCUUU", ".((...))((...))."));

  // Flush coax stack with lonely pairs:
  // Not counted as continuous, so no penultimate stacking.
  EXPECT_EQ(E(13.40), GetEnergy(m, "AGCCCUACCCUU", ".(...)(...)."));

  // Mismatch mediated coax stack
  // GU end on AU at 1, 9; AU end on GU at 2, 10
  EXPECT_EQ(E(9.40 - 0.71 * 2 - 0.31 * 2), GetEnergy(m, "AGACCCUUAGACCCUUU", ".((...)).((...))."));

  // GU end on GC at 2; GU end on AU at 9; AU end on GU at 10
  EXPECT_EQ(E(7.80 + 0.13 - 0.71 - 0.31), GetEnergy(m, "AGGCCCUCAGACCCUUU", ".((...)).((...))."));

  // Mismatch mediated coax stack with lonely pairs:
  // Not counted as continuous, so no penultimate stacking.
  EXPECT_EQ(E(10.60), GetEnergy(m, "AGCCCUAACCCUU", ".(...).(...)."));

  // Other tests:
  EXPECT_EQ(E(4.41), GetEnergy(m, "GCAAAGCC", "((...).)"));
  EXPECT_EQ(E(5.69), GetEnergy(m, "CCCAAAAUG", ".(.(...))"));
  EXPECT_EQ(E(4.50), GetEnergy(m, "UACAGA", "(....)"));
  EXPECT_EQ(E(-1.25), GetEnergy(m, "AGGGUCAUCCG", ".(((...)))."));
  EXPECT_EQ(E(7.00), GetEnergy(m, "AGAGAAACAAAU", "(..(...)...)"));
  EXPECT_EQ(E(9.50), GetEnergy(m, "CGUUGCCUAAAAAGGAAACAAG", "(.............(...)..)"));
  EXPECT_EQ(E(7.70), GetEnergy(m, "CCCGAAACAG", "(..(...).)"));
  EXPECT_EQ(E(7.32), GetEnergy(m, "GACAGAAACGCUGAAUC", "((..(...)......))"));
  EXPECT_EQ(E(16.30), GetEnergy(m, "CUGAAACUGGAAACAGAAAUG", "(.(...)..(...).(...))"));
  EXPECT_EQ(E(17.18), GetEnergy(m, "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))"));
  EXPECT_EQ(E(15.60), GetEnergy(m, "AGCUAAAAACAAAGGUGAAACGU", "(..(...).(...)..(...).)"));
  EXPECT_EQ(E(11.10), GetEnergy(m, "CUGAAACUGGAAACAGAAAUG", ".(.(...)(....)......)"));
  EXPECT_EQ(E(-29.04),
      GetEnergy(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
          "(((((((((((.((...((((....))))..)).)))..((((..((((....))))...)))).))))))))...."));
  EXPECT_EQ(E(14.71), GetEnergy(m, "UCUGAGUAAAUUGCUACGCG", "(....)((...).......)"));
  EXPECT_EQ(E(-45.38), GetEnergy(m, k16sHSapiens3));

  EXPECT_EQ(E(1.61), GetEnergy(m, "GGUCAAAGGUC", "((((...))))"));
  EXPECT_EQ(E(-4.44), GetEnergy(m, "GGGGAAACCCC", "((((...))))"));
  EXPECT_EQ(E(6.20), GetEnergy(m, "UGACAAAGGCGA", "(..(...)...)"));
}

TEST_P(EnergyTestT22, T22P2Pseudofree) {
  auto m = t22_ms[GetParam()];

  {
    // Example from https://doi.org/10.1093/nar/gkac261
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "UGUCGAUACCCUGUCGAUA", "((((((((...))))))))");
    EXPECT_EQ(E(-3.65) + pseudofree, energy);
  }

  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "UAGGUCAGCCCCUGGUCUA", "((((((((...))))))))");
    EXPECT_EQ(E(-4.05) + pseudofree, energy);
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
    EXPECT_EQ(E(4.89 + 0.44) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAGCCCCUU", ".((...)).");
    EXPECT_EQ(E(4.19 + 0.44) + pseudofree, energy);
  }

  // Counting twice, two AU on AUs (0.22):
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AACCCUU", "((...))");
    EXPECT_EQ(E(5.96 + 0.22 * 2) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAACCCUUU", ".((...)).");
    EXPECT_EQ(E(5.26 + 0.22 * 2) + pseudofree, energy);
  }

  // Special hairpin with lonely pair unaffected:
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "ACAGUGCU", "(......)");
    EXPECT_EQ(E(2.40) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAACAGUGCUUU", "..(......)..");
    EXPECT_EQ(E(1.70) + pseudofree, energy);
  }

  // Special hairpins:
  // GU end on AU at 0; AU end on GU at 1
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GACAGUGCUU", "((......))");
    EXPECT_EQ(E(1.80 - 0.31 - 0.71) + pseudofree, energy);
  }
  // GU end on AU at 1; AU end on GU at 2
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGACAGUGCUUU", ".((......)).");
    EXPECT_EQ(E(1.10 - 0.31 - 0.71) + pseudofree, energy);
  }

  // Single nuc bulge loop - treated as continuous.
  // AU end on GU at 0, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGCGAAAAUUUU", "((.((...))))");
    EXPECT_EQ(E(8.42 - 0.71 * 2) + pseudofree, energy);
  }
  // AU end on GU at 1, 5
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GAGCGAAAAUUUUU", ".((.((...)))).");
    EXPECT_EQ(E(7.72 - 0.71 * 2) + pseudofree, energy);
  }

  // Single nuc bulge, continuous single on outer:
  // GU end on GU at 0; AU end on GU at 3
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GCGAAAAUUU", "(.((...)))");
    EXPECT_EQ(E(8.40 - 0.74 - 0.71) + pseudofree, energy);
  }
  // GU end on GU at 1; AU end on GU at 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GGCGAAAAUUUU", ".(.((...))).");
    EXPECT_EQ(E(7.70 - 0.74 - 0.71) + pseudofree, energy);
  }

  // Single nuc bulge, continuous single on inner:
  // AU end on GU at 0, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGCAAAAUUU", "((.(...)))");
    EXPECT_EQ(E(8.62 - 0.71 * 2) + pseudofree, energy);
  }
  // AU end on GU at 1, 5
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GAGCAAAAUUUU", ".((.(...))).");
    EXPECT_EQ(E(7.92 - 0.71 * 2) + pseudofree, energy);
  }

  // Single nuc bulge, continuous single on both:
  // GU end on AU at 0; AU end on GU at 2
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GCAAAAUU", "(.(...))");
    EXPECT_EQ(E(8.60 - 0.31 - 0.71) + pseudofree, energy);
  }
  // GU end on AU at 1; AU end on GU at 3
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GGCAAAAUUU", ".(.(...)).");
    EXPECT_EQ(E(7.90 - 0.31 - 0.71) + pseudofree, energy);
  }

  // Multi nuc bulge loop:
  // AU end on GU at 0, 5; GU end on AU at 1, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGCCGAAAAUUUU", "((..((...))))");
    EXPECT_EQ(E(7.62 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }
  // AU end on GU at 1, 6; GU end on AU at 2, 5
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GAGCCAGAAAUUUUU", ".((..((...)))).");
    EXPECT_EQ(E(7.54 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }

  // Multi nuc bulge loop with lonely pairs:
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GCCGAAAUU", "(..(...))");
    EXPECT_EQ(E(8.20) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GGCCGAAAUUU", ".(..(...)).");
    EXPECT_EQ(E(7.50) + pseudofree, energy);
  }

  // Internal loop:
  // AU end on GU at 0, 4; GU end on AU at 1, 3
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGGGACCCUUCCCUU", "((.((...))...))");
    EXPECT_EQ(E(9.02 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }
  // AU end on GU at 1, 5; GU end on AU at 2, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GAGGGACCCUUCCCUUU", ".((.((...))...)).");
    EXPECT_EQ(E(8.32 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }

  // Internal loop with lonely pairs:
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAACCCUCCCU", "(.(...)...)");
    EXPECT_EQ(E(9.60) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAAACCCUCCCUU", ".(.(...)...).");
    EXPECT_EQ(E(8.90) + pseudofree, energy);
  }

  // Special 1x1 internal loop with helix affected by penultimate stack on external
  // sides only:
  // AU end on GU at 0, 4; GU end on AU at 1, 3
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGGGACCCUUCUU", "((.((...)).))");
    EXPECT_EQ(E(7.22 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }
  // AU end on GU at 1, 5; GU end on AU at 2, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GAGGGACCCUUCUUU", ".((.((...)).)).");
    EXPECT_EQ(E(6.52 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }

  // Special 1x1 internal loop with lonely pairs:
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAACCCUAU", "(.(...).)");
    EXPECT_EQ(E(7.80) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAAACCCUAUU", ".(.(...).).");
    EXPECT_EQ(E(7.10) + pseudofree, energy);
  }

  // Special 1x2 internal loop with helix affected by penultimate stack on external
  // sides only:
  // AU end on GU at 0, 4; GU end on AU at 1, 3
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGGGACCCUUCCUU", "((.((...))..))");
    EXPECT_EQ(E(9.02 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }
  // AU end on GU at 1, 5; GU end on AU at 2, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GAGGGACCCUUCCUUU", ".((.((...))..)).");
    EXPECT_EQ(E(8.32 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }

  // Special 1x2 internal loop with lonely pairs:
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAACCCUCAU", "(.(...)..)");
    EXPECT_EQ(E(9.60) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAAACCCUCAUU", ".(.(...)..).");
    EXPECT_EQ(E(8.90) + pseudofree, energy);
  }

  // Special 2x2 internal loop with helix affected by penultimate stack on external
  // sides only:
  // AU end on GU at 0, 5; GU end on AU at 1, 4
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGGGGACCCUUCCUU", "((..((...))..))");
    EXPECT_EQ(E(8.42 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }
  // AU end on GU at 1, 6; GU end on AU at 2, 5
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GAGGGGACCCUUCCUUU", ".((..((...))..)).");
    EXPECT_EQ(E(7.72 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }

  // Special 2x2 internal loop with lonely pairs:
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAAACCCUAAU", "(..(...)..)");
    EXPECT_EQ(E(8.70) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AAAAACCCUAAUU", ".(..(...)..).");
    EXPECT_EQ(E(8.00) + pseudofree, energy);
  }

  // Multiloop:
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "ACACCCUCCACCCUCCACCCUCU", "(.(...)..(...)..(...).)");
    EXPECT_EQ(E(23.70) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "AACACCCUCCACCCUCCACCCUCUU", ".(.(...)..(...)..(...).).");
    EXPECT_EQ(E(23.00) + pseudofree, energy);
  }

  // AU end on GU at 0; GU end on AU at 1
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "AGCACCCUCCACCCUCCACCCUCUU", "((.(...)..(...)..(...).))");
    EXPECT_EQ(E(23.72 - 0.71 - 0.31) + pseudofree, energy);
  }

  // AU end on GU at 1; GU end on AU at 2
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "AAGCACCCUCCACCCUCCACCCUCUUU", ".((.(...)..(...)..(...).)).");
    EXPECT_EQ(E(23.02 - 0.71 - 0.31) + pseudofree, energy);
  }

  // AU end on GU at 1, 4, 13; GU end on AU at 2, 5, 14
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "AAGCAGCCCUUCCAGCCCUUCCACCCUCUUU", ".((.((...))..((...))..(...).)).");
    EXPECT_EQ(E(23.06 - 0.71 * 3 - 0.31 * 3) + pseudofree, energy);
  }

  // Flush coax stack
  // GU end on AU at 0, 7; AU end on GU at 1, 8 (not taken as continuous)
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GACCCUUGACCCUU", "((...))((...))");
    EXPECT_EQ(E(12.22 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }

  // GU end on GC at 1; GU end on AU at 7; AU end on GU at 8
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GGCCCUCGACCCUU", "((...))((...))");
    EXPECT_EQ(E(10.35 + 0.13 - 0.71 - 0.31) + pseudofree, energy);
  }

  // GU end on AU at 1, 8; AU end on GU at 2, 9
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGACCCUUGACCCUUU", ".((...))((...)).");
    EXPECT_EQ(E(12.20 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }

  // GU end on GC at 2; GU end on AU at 8; AU end on GU at 9
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGGCCCUCGACCCUUU", ".((...))((...)).");
    EXPECT_EQ(E(10.35 + 0.13 - 0.71 - 0.31) + pseudofree, energy);
  }

  // Flush coax stack with lonely pairs:
  // Not counted as continuous, so no penultimate stacking.
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGCCCUACCCUU", ".(...)(...).");
    EXPECT_EQ(E(13.40) + pseudofree, energy);
  }

  // Mismatch mediated coax stack
  // GU end on AU at 1, 9; AU end on GU at 2, 10
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGACCCUUAGACCCUUU", ".((...)).((...)).");
    EXPECT_EQ(E(9.40 - 0.71 * 2 - 0.31 * 2) + pseudofree, energy);
  }

  // GU end on GC at 2; GU end on AU at 9; AU end on GU at 10
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGGCCCUCAGACCCUUU", ".((...)).((...)).");
    EXPECT_EQ(E(7.80 + 0.13 - 0.71 - 0.31) + pseudofree, energy);
  }

  // Mismatch mediated coax stack with lonely pairs:
  // Not counted as continuous, so no penultimate stacking.
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGCCCUAACCCUU", ".(...).(...).");
    EXPECT_EQ(E(10.60) + pseudofree, energy);
  }

  // Other tests:
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GCAAAGCC", "((...).)");
    EXPECT_EQ(E(4.41) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "CCCAAAAUG", ".(.(...))");
    EXPECT_EQ(E(5.69) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "UACAGA", "(....)");
    EXPECT_EQ(E(4.50) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGGGUCAUCCG", ".(((...))).");
    EXPECT_EQ(E(-1.25) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "AGAGAAACAAAU", "(..(...)...)");
    EXPECT_EQ(E(7.00) + pseudofree, energy);
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
    EXPECT_EQ(E(7.32) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "CUGAAACUGGAAACAGAAAUG", "(.(...)..(...).(...))");
    EXPECT_EQ(E(16.30) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))");
    EXPECT_EQ(E(17.18) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "AGCUAAAAACAAAGGUGAAACGU", "(..(...).(...)..(...).)");
    EXPECT_EQ(E(15.60) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "CUGAAACUGGAAACAGAAAUG", ".(.(...)(....)......)");
    EXPECT_EQ(E(11.10) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m,
        "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
        "(((((((((((.((...((((....))))..)).)))..((((..((((....))))...)))).))))))))....");
    EXPECT_EQ(E(-29.04) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] =
        GetPseudofree(m, "UCUGAGUAAAUUGCUACGCG", "(....)((...).......)");
    EXPECT_EQ(E(14.71) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(
        m, std::get<Primary>(k16sHSapiens3).ToSeq(), std::get<Secondary>(k16sHSapiens3).ToDb());
    EXPECT_EQ(E(-45.38) + pseudofree, energy);
  }

  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GGUCAAAGGUC", "((((...))))");
    EXPECT_EQ(E(1.61) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "GGGGAAACCCC", "((((...))))");
    EXPECT_EQ(E(-4.44) + pseudofree, energy);
  }
  {
    const auto& [energy, pseudofree] = GetPseudofree(m, "UGACAAAGGCGA", "(..(...)...)");
    EXPECT_EQ(E(6.20) + pseudofree, energy);
  }
}

#endif

INSTANTIATE_TEST_SUITE_P(EnergyModelTests, EnergyTestT22, testing::Range(0, NUM_T22_MODELS));

}  // namespace mrna
