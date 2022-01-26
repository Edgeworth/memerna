// Copyright 2016 E.
#include "compute/energy/energy.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "common_test.h"
#include "compute/energy/branch.h"
#include "compute/energy/model.h"
#include "compute/energy/precomp.h"
#include "gtest/gtest.h"
#include "model/base.h"
#include "model/model.h"
#include "model/primary.h"
#include "model/secondary.h"

namespace mrna::energy {

class EnergyTest : public testing::Test {
 public:
  std::tuple<Primary, Secondary> kNNDBHairpin1 =
      ParsePrimaryDotBracket("CACAAAAAAAUGUG", "((((......))))");
  std::tuple<Primary, Secondary> kNNDBHairpin2 =
      ParsePrimaryDotBracket("CACAGGAAGUGUG", "((((.....))))");
  std::tuple<Primary, Secondary> kNNDBHairpin3 =
      ParsePrimaryDotBracket("CACCCGAGGGUG", "((((....))))");
  std::tuple<Primary, Secondary> kNNDBHairpin4 =
      ParsePrimaryDotBracket("CACACCCCCCUGUG", "((((......))))");
  std::tuple<Primary, Secondary> kNNDBHairpin5 =
      ParsePrimaryDotBracket("CGGGGGAAGUCCG", "((((.....))))");
  std::tuple<Primary, Secondary> kNNDBBulge1 =
      ParsePrimaryDotBracket("GCCCGAAACGGC", "(((.(...))))");
  std::tuple<Primary, Secondary> kNNDBBulge2 =
      ParsePrimaryDotBracket("GAACAGAAACUC", "((...(...)))");
  std::tuple<Primary, Secondary> kNNDBInternal2x3 =
      ParsePrimaryDotBracket("CAGACGAAACGGAGUG", "((..((...))...))");
  std::tuple<Primary, Secondary> kNNDBInternal1x5 =
      ParsePrimaryDotBracket("CAGCGAAACGGAAAGUG", "((.((...)).....))");
  std::tuple<Primary, Secondary> kNNDBInternal2x2 =
      ParsePrimaryDotBracket("CAGACGAAACGGAUG", "((..((...))..))");
  std::tuple<Primary, Secondary> kFlushCoax =
      ParsePrimaryDotBracket("GUGAAACACAAAAUGA", ".((...))((...)).");
  // NNDB T99 Multiloop example
  std::tuple<Primary, Secondary> kNNDBMultiloop =
      ParsePrimaryDotBracket("UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))");

  std::tuple<Primary, Secondary> kBulge1 = ParsePrimaryDotBracket("GCUCGAAACAGC", "(((.(...))))");
  std::tuple<Primary, Secondary> kInternal1 =
      ParsePrimaryDotBracket("AGAGAAACAAAU", "(..(...)...)");

  Energy GetEnergy(const std::string& r, const std::string& db) {
    return GetEnergy({Primary::FromString(r), Secondary::FromDotBracket(db)});
  }

  Energy GetEnergy(const std::tuple<Primary, Secondary>& s) {
    return g_em.TotalEnergy(std::get<Primary>(s), std::get<Secondary>(s), nullptr).energy;
  }
};

TEST_F(EnergyTest, MultiloopEnergy) {
  EXPECT_EQ(g_em.multiloop_hack_a + 4 * g_em.multiloop_hack_b, g_em.MultiloopInitiation(4));
}

TEST_F(EnergyTest, NNDBHairpinLoopExamples) {
  EXPECT_EQ(g_em.stack[C][A][U][G] + g_em.stack[A][C][G][U] + g_em.stack[C][A][U][G] +
          g_em.augu_penalty + g_em.terminal[A][A][A][U] + g_em.HairpinInitiation(6),
      GetEnergy(kNNDBHairpin1));
  EXPECT_EQ(g_em.stack[C][A][U][G] + g_em.stack[A][C][G][U] + g_em.stack[C][A][U][G] +
          g_em.augu_penalty + g_em.terminal[A][G][G][U] + g_em.hairpin_gg_first_mismatch +
          g_em.HairpinInitiation(5),
      GetEnergy(kNNDBHairpin2));
  EXPECT_EQ(g_em.stack[C][A][U][G] + g_em.stack[A][C][G][U] + g_em.stack[C][C][G][G] +
          g_em.hairpin["CCGAGG"],
      GetEnergy(kNNDBHairpin3));
  EXPECT_EQ(g_em.stack[C][A][U][G] + g_em.stack[A][C][G][U] + g_em.stack[C][A][U][G] +
          g_em.augu_penalty + g_em.terminal[A][C][C][U] + g_em.HairpinInitiation(6) +
          g_em.hairpin_all_c_a * 6 + g_em.hairpin_all_c_b,
      GetEnergy(kNNDBHairpin4));
  EXPECT_EQ(g_em.stack[C][G][C][G] + g_em.stack[G][G][C][C] + g_em.stack[G][G][U][C] +
          g_em.augu_penalty + g_em.terminal[G][G][G][U] + g_em.hairpin_gg_first_mismatch +
          g_em.HairpinInitiation(5) + g_em.hairpin_special_gu_closure,
      GetEnergy(kNNDBHairpin5));

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin1)), g_em);
    EXPECT_EQ(g_em.augu_penalty + g_em.terminal[A][A][A][U] + g_em.HairpinInitiation(6),
        pc.Hairpin(3, 10));
  }

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin2)), g_em);
    EXPECT_EQ(g_em.augu_penalty + g_em.terminal[A][G][G][U] + g_em.hairpin_gg_first_mismatch +
            g_em.HairpinInitiation(5),
        pc.Hairpin(3, 9));
  }

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin3)), g_em);
    EXPECT_EQ(g_em.hairpin["CCGAGG"], pc.Hairpin(3, 8));
  }

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin4)), g_em);
    EXPECT_EQ(g_em.augu_penalty + g_em.terminal[A][C][C][U] + g_em.HairpinInitiation(6) +
            g_em.hairpin_all_c_a * 6 + g_em.hairpin_all_c_b,
        pc.Hairpin(3, 10));
  }

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin5)), g_em);
    EXPECT_EQ(g_em.augu_penalty + g_em.terminal[G][G][G][U] + g_em.hairpin_gg_first_mismatch +
            g_em.HairpinInitiation(5) + g_em.hairpin_special_gu_closure,
        pc.Hairpin(3, 9));
  }
}

TEST_F(EnergyTest, NNDBBulgeLoopExamples) {
  EXPECT_EQ(g_em.stack[G][C][G][C] + g_em.stack[C][C][G][G] + g_em.BulgeInitiation(1) +
          g_em.bulge_special_c + g_em.stack[C][G][C][G] + g_em.HairpinInitiation(3) -
          Energy(round(10.0 * R * T * log(3))),
      GetEnergy(kNNDBBulge1));

  EXPECT_EQ(g_em.stack[G][A][U][C] + g_em.augu_penalty + g_em.BulgeInitiation(3) +
          g_em.HairpinInitiation(3),
      GetEnergy(kNNDBBulge2));
}

TEST_F(EnergyTest, NNDBMultiloopExamples) {
  EXPECT_EQ(g_em.stack[C][A][U][G] + g_em.stack[A][C][G][U] + g_em.stack[C][A][U][G] +
          2 * g_em.augu_penalty + 2 * g_em.HairpinInitiation(3),
      GetEnergy(kFlushCoax));
  EXPECT_EQ(g_em.stack[G][A][U][C] + g_em.terminal[C][G][A][G] + g_em.coax_mismatch_non_contiguous +
          3 * g_em.HairpinInitiation(3) + g_em.MultiloopInitiation(4) + 2 * g_em.augu_penalty,
      GetEnergy(kNNDBMultiloop));
}

TEST_F(EnergyTest, NNDBInternalLoopExamples) {
  EXPECT_EQ(g_em.stack[C][A][U][G] + g_em.stack[C][G][C][G] + g_em.InternalLoopInitiation(5) +
          std::min(g_em.internal_asym, NINIO_MAX_ASYM) + g_em.internal_2x3_mismatch[A][G][G][U] +
          g_em.internal_2x3_mismatch[G][G][A][C] + g_em.internal_augu_penalty +
          g_em.HairpinInitiation(3),
      GetEnergy(kNNDBInternal2x3));
  EXPECT_EQ(g_em.stack[C][A][U][G] + g_em.stack[C][G][C][G] +
          g_em.internal_2x2[A][G][A][C][G][G][A][U] + g_em.HairpinInitiation(3),
      GetEnergy(kNNDBInternal2x2));
  EXPECT_EQ(g_em.stack[C][A][U][G] + g_em.stack[C][G][C][G] + g_em.InternalLoopInitiation(6) +
          std::min(4 * g_em.internal_asym, NINIO_MAX_ASYM) + g_em.internal_augu_penalty +
          g_em.HairpinInitiation(3),
      GetEnergy(kNNDBInternal1x5));
}

TEST_F(EnergyTest, BaseCases) {
  EXPECT_EQ(g_em.augu_penalty + g_em.stack[G][A][U][C] + g_em.hairpin_init[3],
      GetEnergy(ParsePrimaryDotBracket("GAAAAUC", "((...))")));
  EXPECT_EQ(g_em.augu_penalty * 2 + g_em.stack[G][A][U][U] + g_em.hairpin_init[3],
      GetEnergy(ParsePrimaryDotBracket("GAAAAUU", "((...))")));
  EXPECT_EQ(g_em.augu_penalty * 2 + g_em.HairpinInitiation(3) +
          std::min(
              g_em.terminal[U][A][A][A], std::min(g_em.dangle3[U][A][A], g_em.dangle5[U][A][A])),
      GetEnergy(ParsePrimaryDotBracket("AAAAAUA", ".(...).")));
  EXPECT_EQ(g_em.augu_penalty * 2 + g_em.HairpinInitiation(3),
      GetEnergy(ParsePrimaryDotBracket("AAAAU", "(...)")));
  EXPECT_EQ(g_em.stack[G][C][G][C] + g_em.stack[C][U][A][G] + g_em.BulgeInitiation(1) +
          g_em.stack[U][G][C][A] + g_em.HairpinInitiation(3),
      GetEnergy(kBulge1));
  EXPECT_EQ(g_em.InternalLoopInitiation(5) + g_em.internal_asym + g_em.internal_augu_penalty +
          g_em.augu_penalty + g_em.internal_2x3_mismatch[A][G][A][U] +
          g_em.internal_2x3_mismatch[C][A][A][G] + g_em.HairpinInitiation(3),
      GetEnergy(kInternal1));
}

TEST_F(EnergyTest, T04Tests) {
  ONLY_FOR_THIS_MODEL(g_em, T04_MODEL_HASH);

  EXPECT_EQ(88, g_em.HairpinInitiation(87));
  EXPECT_EQ(68, g_em.BulgeInitiation(57));
  EXPECT_EQ(46, g_em.InternalLoopInitiation(67));

  EXPECT_EQ(45, GetEnergy("GCAAAGCC", "((...).)"));
  EXPECT_EQ(57, GetEnergy("CCCAAAAUG", ".(.(...))"));
  EXPECT_EQ(55, GetEnergy("UACAGA", "(....)"));
  EXPECT_EQ(-6, GetEnergy("AGGGUCAUCCG", ".(((...)))."));
  EXPECT_EQ(80, GetEnergy("AGAGAAACAAAU", "(..(...)...)"));
  EXPECT_EQ(95, GetEnergy("CGUUGCCUAAAAAGGAAACAAG", "(.............(...)..)"));
  EXPECT_EQ(77, GetEnergy("CCCGAAACAG", "(..(...).)"));
  EXPECT_EQ(74, GetEnergy("GACAGAAACGCUGAAUC", "((..(...)......))"));
  EXPECT_EQ(173, GetEnergy("CUGAAACUGGAAACAGAAAUG", "(.(...)..(...).(...))"));
  EXPECT_EQ(182, GetEnergy("UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))"));
  EXPECT_EQ(176, GetEnergy("AGCUAAAAACAAAGGUGAAACGU", "(..(...).(...)..(...).)"));
  EXPECT_EQ(131, GetEnergy("CUGAAACUGGAAACAGAAAUG", ".(.(...)(....)......)"));
  EXPECT_EQ(-276,
      GetEnergy("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
          "(((((((((((.((...((((....))))..)).)))..((((..((((....))))...)))).))))))))...."));
  EXPECT_EQ(179, GetEnergy("UCUGAGUAAAUUGCUACGCG", "(....)((...).......)"));

  // Special stacking - this is not implemented. TODO: Implement this?
  EXPECT_EQ(37, GetEnergy("GGUCAAAGGUC", "((((...))))"));
  EXPECT_EQ(-45, GetEnergy("GGGGAAACCCC", "((((...))))"));
  EXPECT_EQ(72, GetEnergy("UGACAAAGGCGA", "(..(...)...)"));
}

TEST_F(EnergyTest, Precomp) {
  ONLY_FOR_THIS_MODEL(g_em, T04_MODEL_HASH);

  const Precomp pc(Primary::FromString("GGGGAAACCCC"), g_em);
  EXPECT_EQ(-21 - 4 - 16, pc.min_mismatch_coax);
  EXPECT_EQ(-34, pc.min_flush_coax);
  EXPECT_EQ(-26, pc.min_twoloop_not_stack);

  Energy augubranch[4][4] = {
      {-6, -6, -6, 5 - 6}, {-6, -6, -6, -6}, {-6, -6, -6, 5 - 6}, {5 - 6, -6, 5 - 6, -6}};
  EXPECT_EQ(sizeof(augubranch), sizeof(pc.augubranch));
  EXPECT_EQ(0, std::memcmp(augubranch, pc.augubranch, sizeof(augubranch)));
}

TEST_F(EnergyTest, Helpers) {
  EXPECT_EQ(0, internal::MaxNumContiguous(Primary::FromString("")));
  EXPECT_EQ(1, internal::MaxNumContiguous(Primary::FromString("A")));
  EXPECT_EQ(2, internal::MaxNumContiguous(Primary::FromString("AA")));
  EXPECT_EQ(2, internal::MaxNumContiguous(Primary::FromString("GUAAC")));
  EXPECT_EQ(1, internal::MaxNumContiguous(Primary::FromString("GUACA")));
  EXPECT_EQ(3, internal::MaxNumContiguous(Primary::FromString("GAUCCC")));
  EXPECT_EQ(3, internal::MaxNumContiguous(Primary::FromString("GGGAUC")));
  EXPECT_EQ(4, internal::MaxNumContiguous(Primary::FromString("GGGAUCAAAA")));
  EXPECT_EQ(5, internal::MaxNumContiguous(Primary::FromString("GGGAUUUUUCAAAA")));
}

TEST_F(EnergyTest, GetBranchCounts) {
  EXPECT_EQ((std::vector<int>{2, 0}), GetBranchCounts(Secondary::FromDotBracket("()")));
  EXPECT_EQ((std::vector<int>{}), GetBranchCounts(Secondary::FromDotBracket("")));
  EXPECT_EQ((std::vector<int>{0}), GetBranchCounts(Secondary::FromDotBracket(".")));
  EXPECT_EQ(
      (std::vector<int>{2, 0, 2, 0, 2, 0}), GetBranchCounts(Secondary::FromDotBracket("()()()")));
  EXPECT_EQ((std::vector<int>{2, 1, 0, 1}), GetBranchCounts(Secondary::FromDotBracket("(())")));
  EXPECT_EQ((std::vector<int>{2, 3, 0, 3, 0, 3, 0, 3}),
      GetBranchCounts(Secondary::FromDotBracket("(()()())")));
  EXPECT_EQ((std::vector<int>{2, 1, 2, 0, 2, 0, 2, 1}),
      GetBranchCounts(Secondary::FromDotBracket("((()()))")));
}

}  // namespace mrna::energy
