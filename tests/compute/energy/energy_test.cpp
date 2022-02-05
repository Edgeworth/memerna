// Copyright 2016 Eliot Courtney.
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
      ParseSeqDb("CACAAAAAAAUGUG", "((((......))))");
  std::tuple<Primary, Secondary> kNNDBHairpin2 =
      ParseSeqDb("CACAGGAAGUGUG", "((((.....))))");
  std::tuple<Primary, Secondary> kNNDBHairpin3 =
      ParseSeqDb("CACCCGAGGGUG", "((((....))))");
  std::tuple<Primary, Secondary> kNNDBHairpin4 =
      ParseSeqDb("CACACCCCCCUGUG", "((((......))))");
  std::tuple<Primary, Secondary> kNNDBHairpin5 =
      ParseSeqDb("CGGGGGAAGUCCG", "((((.....))))");
  std::tuple<Primary, Secondary> kNNDBBulge1 =
      ParseSeqDb("GCCCGAAACGGC", "(((.(...))))");
  std::tuple<Primary, Secondary> kNNDBBulge2 =
      ParseSeqDb("GAACAGAAACUC", "((...(...)))");
  std::tuple<Primary, Secondary> kNNDBInternal2x3 =
      ParseSeqDb("CAGACGAAACGGAGUG", "((..((...))...))");
  std::tuple<Primary, Secondary> kNNDBInternal1x5 =
      ParseSeqDb("CAGCGAAACGGAAAGUG", "((.((...)).....))");
  std::tuple<Primary, Secondary> kNNDBInternal2x2 =
      ParseSeqDb("CAGACGAAACGGAUG", "((..((...))..))");
  std::tuple<Primary, Secondary> kFlushCoax =
      ParseSeqDb("GUGAAACACAAAAUGA", ".((...))((...)).");
  // NNDB T99 Multiloop example
  std::tuple<Primary, Secondary> kNNDBMultiloop =
      ParseSeqDb("UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))");

  std::tuple<Primary, Secondary> kBulge1 = ParseSeqDb("GCUCGAAACAGC", "(((.(...))))");
  std::tuple<Primary, Secondary> kInternal1 =
      ParseSeqDb("AGAGAAACAAAU", "(..(...)...)");

  Energy GetEnergy(const std::string& r, const std::string& db) {
    return GetEnergy({Primary::FromSeq(r), Secondary::FromDb(db)});
  }

  Energy GetEnergy(const std::tuple<Primary, Secondary>& s) {
    return t04->TotalEnergy(std::get<Primary>(s), std::get<Secondary>(s), nullptr).energy;
  }
};

TEST_F(EnergyTest, MultiloopEnergy) {
  EXPECT_EQ(t04->multiloop_hack_a + 4 * t04->multiloop_hack_b, t04->MultiloopInitiation(4));
}

TEST_F(EnergyTest, NNDBHairpinLoopExamples) {
  EXPECT_EQ(t04->stack[C][A][U][G] + t04->stack[A][C][G][U] + t04->stack[C][A][U][G] +
          t04->augu_penalty + t04->terminal[A][A][A][U] + t04->HairpinInitiation(6),
      GetEnergy(kNNDBHairpin1));
  EXPECT_EQ(t04->stack[C][A][U][G] + t04->stack[A][C][G][U] + t04->stack[C][A][U][G] +
          t04->augu_penalty + t04->terminal[A][G][G][U] + t04->hairpin_gg_first_mismatch +
          t04->HairpinInitiation(5),
      GetEnergy(kNNDBHairpin2));
  EXPECT_EQ(t04->stack[C][A][U][G] + t04->stack[A][C][G][U] + t04->stack[C][C][G][G] +
          t04->hairpin["CCGAGG"],
      GetEnergy(kNNDBHairpin3));
  EXPECT_EQ(t04->stack[C][A][U][G] + t04->stack[A][C][G][U] + t04->stack[C][A][U][G] +
          t04->augu_penalty + t04->terminal[A][C][C][U] + t04->HairpinInitiation(6) +
          t04->hairpin_all_c_a * 6 + t04->hairpin_all_c_b,
      GetEnergy(kNNDBHairpin4));
  EXPECT_EQ(t04->stack[C][G][C][G] + t04->stack[G][G][C][C] + t04->stack[G][G][U][C] +
          t04->augu_penalty + t04->terminal[G][G][G][U] + t04->hairpin_gg_first_mismatch +
          t04->HairpinInitiation(5) + t04->hairpin_special_gu_closure,
      GetEnergy(kNNDBHairpin5));

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin1)), t04);
    EXPECT_EQ(t04->augu_penalty + t04->terminal[A][A][A][U] + t04->HairpinInitiation(6),
        pc.Hairpin(3, 10));
  }

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin2)), t04);
    EXPECT_EQ(t04->augu_penalty + t04->terminal[A][G][G][U] + t04->hairpin_gg_first_mismatch +
            t04->HairpinInitiation(5),
        pc.Hairpin(3, 9));
  }

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin3)), t04);
    EXPECT_EQ(t04->hairpin["CCGAGG"], pc.Hairpin(3, 8));
  }

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin4)), t04);
    EXPECT_EQ(t04->augu_penalty + t04->terminal[A][C][C][U] + t04->HairpinInitiation(6) +
            t04->hairpin_all_c_a * 6 + t04->hairpin_all_c_b,
        pc.Hairpin(3, 10));
  }

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin5)), t04);
    EXPECT_EQ(t04->augu_penalty + t04->terminal[G][G][G][U] + t04->hairpin_gg_first_mismatch +
            t04->HairpinInitiation(5) + t04->hairpin_special_gu_closure,
        pc.Hairpin(3, 9));
  }
}

TEST_F(EnergyTest, NNDBBulgeLoopExamples) {
  EXPECT_EQ(t04->stack[G][C][G][C] + t04->stack[C][C][G][G] + t04->BulgeInitiation(1) +
          t04->bulge_special_c + t04->stack[C][G][C][G] + t04->HairpinInitiation(3) -
          Energy(round(10.0 * R * T * log(3))),
      GetEnergy(kNNDBBulge1));

  EXPECT_EQ(t04->stack[G][A][U][C] + t04->augu_penalty + t04->BulgeInitiation(3) +
          t04->HairpinInitiation(3),
      GetEnergy(kNNDBBulge2));
}

TEST_F(EnergyTest, NNDBMultiloopExamples) {
  EXPECT_EQ(t04->stack[C][A][U][G] + t04->stack[A][C][G][U] + t04->stack[C][A][U][G] +
          2 * t04->augu_penalty + 2 * t04->HairpinInitiation(3),
      GetEnergy(kFlushCoax));
  EXPECT_EQ(t04->stack[G][A][U][C] + t04->terminal[C][G][A][G] + t04->coax_mismatch_non_contiguous +
          3 * t04->HairpinInitiation(3) + t04->MultiloopInitiation(4) + 2 * t04->augu_penalty,
      GetEnergy(kNNDBMultiloop));
}

TEST_F(EnergyTest, NNDBInternalLoopExamples) {
  EXPECT_EQ(t04->stack[C][A][U][G] + t04->stack[C][G][C][G] + t04->InternalLoopInitiation(5) +
          std::min(t04->internal_asym, NINIO_MAX_ASYM) + t04->internal_2x3_mismatch[A][G][G][U] +
          t04->internal_2x3_mismatch[G][G][A][C] + t04->internal_augu_penalty +
          t04->HairpinInitiation(3),
      GetEnergy(kNNDBInternal2x3));
  EXPECT_EQ(t04->stack[C][A][U][G] + t04->stack[C][G][C][G] +
          t04->internal_2x2[A][G][A][C][G][G][A][U] + t04->HairpinInitiation(3),
      GetEnergy(kNNDBInternal2x2));
  EXPECT_EQ(t04->stack[C][A][U][G] + t04->stack[C][G][C][G] + t04->InternalLoopInitiation(6) +
          std::min(4 * t04->internal_asym, NINIO_MAX_ASYM) + t04->internal_augu_penalty +
          t04->HairpinInitiation(3),
      GetEnergy(kNNDBInternal1x5));
}

TEST_F(EnergyTest, BaseCases) {
  EXPECT_EQ(t04->augu_penalty + t04->stack[G][A][U][C] + t04->hairpin_init[3],
      GetEnergy(ParseSeqDb("GAAAAUC", "((...))")));
  EXPECT_EQ(t04->augu_penalty * 2 + t04->stack[G][A][U][U] + t04->hairpin_init[3],
      GetEnergy(ParseSeqDb("GAAAAUU", "((...))")));
  EXPECT_EQ(t04->augu_penalty * 2 + t04->HairpinInitiation(3) +
          std::min(
              t04->terminal[U][A][A][A], std::min(t04->dangle3[U][A][A], t04->dangle5[U][A][A])),
      GetEnergy(ParseSeqDb("AAAAAUA", ".(...).")));
  EXPECT_EQ(t04->augu_penalty * 2 + t04->HairpinInitiation(3),
      GetEnergy(ParseSeqDb("AAAAU", "(...)")));
  EXPECT_EQ(t04->stack[G][C][G][C] + t04->stack[C][U][A][G] + t04->BulgeInitiation(1) +
          t04->stack[U][G][C][A] + t04->HairpinInitiation(3),
      GetEnergy(kBulge1));
  EXPECT_EQ(t04->InternalLoopInitiation(5) + t04->internal_asym + t04->internal_augu_penalty +
          t04->augu_penalty + t04->internal_2x3_mismatch[A][G][A][U] +
          t04->internal_2x3_mismatch[C][A][A][G] + t04->HairpinInitiation(3),
      GetEnergy(kInternal1));
}

TEST_F(EnergyTest, T04Tests) {
  EXPECT_EQ(88, t04->HairpinInitiation(87));
  EXPECT_EQ(68, t04->BulgeInitiation(57));
  EXPECT_EQ(46, t04->InternalLoopInitiation(67));

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
  const Precomp pc(Primary::FromSeq("GGGGAAACCCC"), t04);
  EXPECT_EQ(-21 - 4 - 16, pc.min_mismatch_coax);
  EXPECT_EQ(-34, pc.min_flush_coax);
  EXPECT_EQ(-26, pc.min_twoloop_not_stack);

  Energy augubranch[4][4] = {
      {-6, -6, -6, 5 - 6}, {-6, -6, -6, -6}, {-6, -6, -6, 5 - 6}, {5 - 6, -6, 5 - 6, -6}};
  EXPECT_EQ(sizeof(augubranch), sizeof(pc.augubranch));
  EXPECT_EQ(0, std::memcmp(augubranch, pc.augubranch, sizeof(augubranch)));
}

TEST_F(EnergyTest, Helpers) {
  EXPECT_EQ(0, internal::MaxNumContiguous(Primary::FromSeq("")));
  EXPECT_EQ(1, internal::MaxNumContiguous(Primary::FromSeq("A")));
  EXPECT_EQ(2, internal::MaxNumContiguous(Primary::FromSeq("AA")));
  EXPECT_EQ(2, internal::MaxNumContiguous(Primary::FromSeq("GUAAC")));
  EXPECT_EQ(1, internal::MaxNumContiguous(Primary::FromSeq("GUACA")));
  EXPECT_EQ(3, internal::MaxNumContiguous(Primary::FromSeq("GAUCCC")));
  EXPECT_EQ(3, internal::MaxNumContiguous(Primary::FromSeq("GGGAUC")));
  EXPECT_EQ(4, internal::MaxNumContiguous(Primary::FromSeq("GGGAUCAAAA")));
  EXPECT_EQ(5, internal::MaxNumContiguous(Primary::FromSeq("GGGAUUUUUCAAAA")));
}

TEST_F(EnergyTest, GetBranchCounts) {
  EXPECT_EQ((std::vector<int>{2, 0}), GetBranchCounts(Secondary::FromDb("()")));
  EXPECT_EQ((std::vector<int>{}), GetBranchCounts(Secondary::FromDb("")));
  EXPECT_EQ((std::vector<int>{0}), GetBranchCounts(Secondary::FromDb(".")));
  EXPECT_EQ(
      (std::vector<int>{2, 0, 2, 0, 2, 0}), GetBranchCounts(Secondary::FromDb("()()()")));
  EXPECT_EQ((std::vector<int>{2, 1, 0, 1}), GetBranchCounts(Secondary::FromDb("(())")));
  EXPECT_EQ((std::vector<int>{2, 3, 0, 3, 0, 3, 0, 3}),
      GetBranchCounts(Secondary::FromDb("(()()())")));
  EXPECT_EQ((std::vector<int>{2, 1, 2, 0, 2, 0, 2, 1}),
      GetBranchCounts(Secondary::FromDb("((()()))")));
}

}  // namespace mrna::energy
