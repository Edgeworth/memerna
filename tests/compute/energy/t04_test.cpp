// Copyright 2016 Eliot Courtney.
#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "common_test.h"
#include "compute/energy/branch.h"
#include "compute/energy/energy.h"
#include "compute/energy/t04/precomp.h"
#include "gtest/gtest.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/primary.h"
#include "model/secondary.h"

namespace mrna::energy {

class T04LikeModelTest : public testing::TestWithParam<t04::ModelPtr> {
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
  std::tuple<Primary, Secondary> kFlushCoax = ParseSeqDb("GUGAAACACAAAAUGA", ".((...))((...)).");
  // NNDB T99 Multiloop example
  std::tuple<Primary, Secondary> kNNDBMultiloop =
      ParseSeqDb("UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))");

  std::tuple<Primary, Secondary> kBulge1 = ParseSeqDb("GCUCGAAACAGC", "(((.(...))))");
  std::tuple<Primary, Secondary> kInternal1 = ParseSeqDb("AGAGAAACAAAU", "(..(...)...)");

  static Energy GetEnergy(const std::string& r, const std::string& db) {
    return GetEnergy({Primary::FromSeq(r), Secondary::FromDb(db)});
  }

  static Energy GetEnergy(const std::tuple<Primary, Secondary>& s) {
    return GetParam()->TotalEnergy(std::get<Primary>(s), std::get<Secondary>(s), nullptr).energy;
  }
};

TEST_P(T04LikeModelTest, MultiloopEnergy) {
  auto em = GetParam();
  EXPECT_EQ(em->multiloop_hack_a + 4 * em->multiloop_hack_b, em->MultiloopInitiation(4));
}

TEST_P(T04LikeModelTest, NNDBHairpinLoopExamples) {
  auto em = GetParam();

  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][A][U][G] +
          em->augu_penalty + em->terminal[A][A][A][U] + em->HairpinInitiation(6),
      GetEnergy(kNNDBHairpin1));
  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][A][U][G] +
          em->augu_penalty + em->terminal[A][G][G][U] + em->hairpin_gg_first_mismatch +
          em->HairpinInitiation(5),
      GetEnergy(kNNDBHairpin2));
  EXPECT_EQ(
      em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][C][G][G] + em->hairpin["CCGAGG"],
      GetEnergy(kNNDBHairpin3));
  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][A][U][G] +
          em->augu_penalty + em->terminal[A][C][C][U] + em->HairpinInitiation(6) +
          em->hairpin_all_c_a * 6 + em->hairpin_all_c_b,
      GetEnergy(kNNDBHairpin4));
  EXPECT_EQ(em->stack[C][G][C][G] + em->stack[G][G][C][C] + em->stack[G][G][U][C] +
          em->augu_penalty + em->terminal[G][G][G][U] + em->hairpin_gg_first_mismatch +
          em->HairpinInitiation(5) + em->hairpin_special_gu_closure,
      GetEnergy(kNNDBHairpin5));

  {
    const t04::Precomp pc(Primary(std::get<Primary>(kNNDBHairpin1)), em);
    EXPECT_EQ(
        em->augu_penalty + em->terminal[A][A][A][U] + em->HairpinInitiation(6), pc.Hairpin(3, 10));
  }

  {
    const t04::Precomp pc(Primary(std::get<Primary>(kNNDBHairpin2)), em);
    EXPECT_EQ(em->augu_penalty + em->terminal[A][G][G][U] + em->hairpin_gg_first_mismatch +
            em->HairpinInitiation(5),
        pc.Hairpin(3, 9));
  }

  {
    const t04::Precomp pc(Primary(std::get<Primary>(kNNDBHairpin3)), em);
    EXPECT_EQ(em->hairpin["CCGAGG"], pc.Hairpin(3, 8));
  }

  {
    const t04::Precomp pc(Primary(std::get<Primary>(kNNDBHairpin4)), em);
    EXPECT_EQ(em->augu_penalty + em->terminal[A][C][C][U] + em->HairpinInitiation(6) +
            em->hairpin_all_c_a * 6 + em->hairpin_all_c_b,
        pc.Hairpin(3, 10));
  }

  {
    const t04::Precomp pc(Primary(std::get<Primary>(kNNDBHairpin5)), em);
    EXPECT_EQ(em->augu_penalty + em->terminal[G][G][G][U] + em->hairpin_gg_first_mismatch +
            em->HairpinInitiation(5) + em->hairpin_special_gu_closure,
        pc.Hairpin(3, 9));
  }
}

TEST_P(T04LikeModelTest, NNDBBulgeLoopExamples) {
  auto em = GetParam();

  EXPECT_EQ(em->stack[G][C][G][C] + em->stack[C][C][G][G] + em->BulgeInitiation(1) +
          em->bulge_special_c + em->stack[C][G][C][G] + em->HairpinInitiation(3) -
          E(R * T * log(3)),
      GetEnergy(kNNDBBulge1));

  EXPECT_EQ(
      em->stack[G][A][U][C] + em->augu_penalty + em->BulgeInitiation(3) + em->HairpinInitiation(3),
      GetEnergy(kNNDBBulge2));
}

TEST_P(T04LikeModelTest, NNDBMultiloopExamples) {
  auto em = GetParam();

  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][A][U][G] +
          2 * em->augu_penalty + 2 * em->HairpinInitiation(3),
      GetEnergy(kFlushCoax));
  EXPECT_EQ(em->stack[G][A][U][C] + em->terminal[C][G][A][G] + em->coax_mismatch_non_contiguous +
          3 * em->HairpinInitiation(3) + em->MultiloopInitiation(4) + 2 * em->augu_penalty,
      GetEnergy(kNNDBMultiloop));
}

TEST_P(T04LikeModelTest, NNDBInternalLoopExamples) {
  auto em = GetParam();

  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[C][G][C][G] + em->InternalLoopInitiation(5) +
          std::min(em->internal_asym, NINIO_MAX_ASYM) + em->internal_2x3_mismatch[A][G][G][U] +
          em->internal_2x3_mismatch[G][G][A][C] + em->internal_augu_penalty +
          em->HairpinInitiation(3),
      GetEnergy(kNNDBInternal2x3));
  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[C][G][C][G] +
          em->internal_2x2[A][G][A][C][G][G][A][U] + em->HairpinInitiation(3),
      GetEnergy(kNNDBInternal2x2));
  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[C][G][C][G] + em->InternalLoopInitiation(6) +
          std::min(4 * em->internal_asym, NINIO_MAX_ASYM) + em->internal_augu_penalty +
          em->HairpinInitiation(3),
      GetEnergy(kNNDBInternal1x5));
}

TEST_P(T04LikeModelTest, BaseCases) {
  auto em = GetParam();

  EXPECT_EQ(em->augu_penalty + em->stack[G][A][U][C] + em->hairpin_init[3],
      GetEnergy("GAAAAUC", "((...))"));
  EXPECT_EQ(em->augu_penalty * 2 + em->stack[G][A][U][U] + em->hairpin_init[3],
      GetEnergy("GAAAAUU", "((...))"));
  EXPECT_EQ(em->augu_penalty * 2 + em->HairpinInitiation(3) +
          std::min(em->terminal[U][A][A][A], std::min(em->dangle3[U][A][A], em->dangle5[U][A][A])),
      GetEnergy("AAAAAUA", ".(...)."));
  EXPECT_EQ(em->augu_penalty * 2 + em->HairpinInitiation(3), GetEnergy("AAAAU", "(...)"));
  EXPECT_EQ(em->stack[G][C][G][C] + em->stack[C][U][A][G] + em->BulgeInitiation(1) +
          em->stack[U][G][C][A] + em->HairpinInitiation(3),
      GetEnergy(kBulge1));
  EXPECT_EQ(em->InternalLoopInitiation(5) + em->internal_asym + em->internal_augu_penalty +
          em->augu_penalty + em->internal_2x3_mismatch[A][G][A][U] +
          em->internal_2x3_mismatch[C][A][A][G] + em->HairpinInitiation(3),
      GetEnergy(kInternal1));
}

TEST_P(T04LikeModelTest, T04Tests) {
  auto em = GetParam();

  EXPECT_EQ(E(8.8), em->HairpinInitiation(87));
  EXPECT_EQ(E(6.8), em->BulgeInitiation(57));
  EXPECT_EQ(E(4.6), em->InternalLoopInitiation(67));

  EXPECT_EQ(E(4.5), GetEnergy("GCAAAGCC", "((...).)"));
  EXPECT_EQ(E(5.7), GetEnergy("CCCAAAAUG", ".(.(...))"));
  EXPECT_EQ(E(5.5), GetEnergy("UACAGA", "(....)"));
  EXPECT_EQ(E(-0.6), GetEnergy("AGGGUCAUCCG", ".(((...)))."));
  EXPECT_EQ(E(8.0), GetEnergy("AGAGAAACAAAU", "(..(...)...)"));
  EXPECT_EQ(E(9.5), GetEnergy("CGUUGCCUAAAAAGGAAACAAG", "(.............(...)..)"));
  EXPECT_EQ(E(7.7), GetEnergy("CCCGAAACAG", "(..(...).)"));
  EXPECT_EQ(E(7.4), GetEnergy("GACAGAAACGCUGAAUC", "((..(...)......))"));
  EXPECT_EQ(E(17.3), GetEnergy("CUGAAACUGGAAACAGAAAUG", "(.(...)..(...).(...))"));
  EXPECT_EQ(E(18.2), GetEnergy("UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))"));
  EXPECT_EQ(E(17.6), GetEnergy("AGCUAAAAACAAAGGUGAAACGU", "(..(...).(...)..(...).)"));
  EXPECT_EQ(E(13.1), GetEnergy("CUGAAACUGGAAACAGAAAUG", ".(.(...)(....)......)"));
  EXPECT_EQ(E(-27.6),
      GetEnergy("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
          "(((((((((((.((...((((....))))..)).)))..((((..((((....))))...)))).))))))))...."));
  EXPECT_EQ(E(17.9), GetEnergy("UCUGAGUAAAUUGCUACGCG", "(....)((...).......)"));

  // Special stacking - this is not implemented. TODO(4): Implement this?
  EXPECT_EQ(E(3.7), GetEnergy("GGUCAAAGGUC", "((((...))))"));
  EXPECT_EQ(E(-4.5), GetEnergy("GGGGAAACCCC", "((((...))))"));
  EXPECT_EQ(E(7.2), GetEnergy("UGACAAAGGCGA", "(..(...)...)"));
}

TEST_P(T04LikeModelTest, Precomp) {
  auto em = GetParam();

  const t04::Precomp pc(Primary::FromSeq("GGGGAAACCCC"), em);
  EXPECT_EQ(E(-2.1 - 0.4 - 1.6), pc.min_mismatch_coax);
  EXPECT_EQ(E(-3.4), pc.min_flush_coax);
  EXPECT_EQ(E(-2.6), pc.min_twoloop_not_stack);

  Energy augubranch[4][4] = {{E(-0.6), E(-0.6), E(-0.6), E(0.5 - 0.6)},
      {E(-0.6), E(-0.6), E(-0.6), E(-0.6)}, {E(-0.6), E(-0.6), E(-0.6), E(0.5 - 0.6)},
      {E(0.5 - 0.6), E(-0.6), E(0.5 - 0.6), E(-0.6)}};
  EXPECT_EQ(sizeof(augubranch), sizeof(pc.augubranch));
  EXPECT_EQ(0, std::memcmp(augubranch, pc.augubranch, sizeof(augubranch)));
}

}  // namespace mrna::energy
