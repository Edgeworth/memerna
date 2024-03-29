// Copyright 2016 E.
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

#include "api/energy/energy.h"
#include "common_test.h"
#include "gtest/gtest.h"
#include "model/base.h"
#include "model/constants.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "models/common/branch.h"
#include "models/t04/energy/branch.h"
#include "models/t04/energy/model.h"
#include "models/t04/energy/precomp.h"

namespace mrna::md::t04 {

Energy GetEnergy(const Model::Ptr& em, const std::tuple<Primary, Secondary>& s) {
  return em->TotalEnergy(std::get<Primary>(s), std::get<Secondary>(s), nullptr).energy;
}

Energy GetEnergy(const Model::Ptr& em, const std::string& r, const std::string& db) {
  return GetEnergy(em, {Primary::FromSeq(r), Secondary::FromDb(db)});
}

class T04ModelTest : public testing::TestWithParam<int> {
 public:
  static Energy GetEnergy(const std::string& r, const std::string& db) {
    return GetEnergy({Primary::FromSeq(r), Secondary::FromDb(db)});
  }

  static Energy GetEnergy(const std::tuple<Primary, Secondary>& s) {
    return test_t04_ems[GetParam()]
        ->TotalEnergy(std::get<Primary>(s), std::get<Secondary>(s), nullptr)
        .energy;
  }
};

TEST_P(T04ModelTest, MultiloopEnergy) {
  auto em = test_t04_ems[GetParam()];
  EXPECT_EQ(em->multiloop_hack_a + 4 * em->multiloop_hack_b, em->MultiloopInitiation(4));
}

TEST_P(T04ModelTest, NNDBHairpinLoopExamples) {
  auto em = test_t04_ems[GetParam()];

  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][A][U][G] + em->au_penalty +
          em->terminal[A][A][A][U] + em->HairpinInitiation(6),
      GetEnergy(kNNDBHairpin1));
  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][A][U][G] + em->au_penalty +
          em->terminal[A][G][G][U] + em->hairpin_gg_first_mismatch + em->HairpinInitiation(5),
      GetEnergy(kNNDBHairpin2));

  if (em->hairpin.contains("CCGAGG")) {
    EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][C][G][G] +
            em->hairpin["CCGAGG"],
        GetEnergy(kNNDBHairpin3));
  }

  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][A][U][G] + em->au_penalty +
          em->terminal[A][C][C][U] + em->HairpinInitiation(6) + em->hairpin_all_c_a * 6 +
          em->hairpin_all_c_b,
      GetEnergy(kNNDBHairpin4));
  EXPECT_EQ(em->stack[C][G][C][G] + em->stack[G][G][C][C] + em->stack[G][G][U][C] + em->gu_penalty +
          em->terminal[G][G][G][U] + em->hairpin_gg_first_mismatch + em->HairpinInitiation(5) +
          em->hairpin_special_gu_closure,
      GetEnergy(kNNDBHairpin5));

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin1)), em);
    EXPECT_EQ(
        em->au_penalty + em->terminal[A][A][A][U] + em->HairpinInitiation(6), pc.Hairpin(3, 10));
  }

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin2)), em);
    EXPECT_EQ(em->au_penalty + em->terminal[A][G][G][U] + em->hairpin_gg_first_mismatch +
            em->HairpinInitiation(5),
        pc.Hairpin(3, 9));
  }

  if (em->hairpin.contains("CCGAGG")) {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin3)), em);
    EXPECT_EQ(em->hairpin["CCGAGG"], pc.Hairpin(3, 8));
  }

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin4)), em);
    EXPECT_EQ(em->au_penalty + em->terminal[A][C][C][U] + em->HairpinInitiation(6) +
            em->hairpin_all_c_a * 6 + em->hairpin_all_c_b,
        pc.Hairpin(3, 10));
  }

  {
    const Precomp pc(Primary(std::get<Primary>(kNNDBHairpin5)), em);
    EXPECT_EQ(em->gu_penalty + em->terminal[G][G][G][U] + em->hairpin_gg_first_mismatch +
            em->HairpinInitiation(5) + em->hairpin_special_gu_closure,
        pc.Hairpin(3, 9));
  }
}

TEST_P(T04ModelTest, NNDBBulgeLoopExamples) {
  auto em = test_t04_ems[GetParam()];

  EXPECT_EQ(em->stack[G][C][G][C] + em->stack[C][C][G][G] + em->BulgeInitiation(1) +
          em->bulge_special_c + em->stack[C][G][C][G] + em->HairpinInitiation(3) -
          E(R * T * log(3)),
      GetEnergy(kNNDBBulge1));

  EXPECT_EQ(
      em->stack[G][A][U][C] + em->au_penalty + em->BulgeInitiation(3) + em->HairpinInitiation(3),
      GetEnergy(kNNDBBulge2));
}

TEST_P(T04ModelTest, NNDBInternalLoopExamples) {
  auto em = test_t04_ems[GetParam()];

  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[C][G][C][G] + em->InternalLoopInitiation(5) +
          std::min(em->internal_asym, NINIO_MAX_ASYM) + em->internal_2x3_mismatch[A][G][G][U] +
          em->internal_2x3_mismatch[G][G][A][C] + em->au_penalty + em->internal_au_penalty +
          em->HairpinInitiation(3),
      GetEnergy(kNNDBInternal2x3));
  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[C][G][C][G] +
          em->internal_2x2[A][G][A][C][G][G][A][U] + em->au_penalty + em->HairpinInitiation(3),
      GetEnergy(kNNDBInternal2x2));
  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[C][G][C][G] + em->InternalLoopInitiation(6) +
          std::min(4 * em->internal_asym, NINIO_MAX_ASYM) + em->au_penalty +
          em->internal_au_penalty + em->HairpinInitiation(3),
      GetEnergy(kNNDBInternal1x5));
}

TEST_P(T04ModelTest, BaseCases) {
  auto em = test_t04_ems[GetParam()];

  EXPECT_EQ(em->au_penalty + em->stack[G][A][U][C] + em->hairpin_init[3],
      GetEnergy("GAAAAUC", "((...))"));
  EXPECT_EQ(em->au_penalty + em->gu_penalty + em->stack[G][A][U][U] + em->hairpin_init[3],
      GetEnergy("GAAAAUU", "((...))"));
  EXPECT_EQ(em->au_penalty * 2 + em->HairpinInitiation(3) +
          std::min(ZERO_E,
              std::min(
                  em->terminal[U][A][A][A], std::min(em->dangle3[U][A][A], em->dangle5[U][A][A]))),
      GetEnergy("AAAAAUA", ".(...)."));
  EXPECT_EQ(em->au_penalty * 2 + em->HairpinInitiation(3), GetEnergy("AAAAU", "(...)"));
  EXPECT_EQ(em->stack[G][C][G][C] + em->stack[C][U][A][G] + em->BulgeInitiation(1) +
          em->stack[U][G][C][A] + em->HairpinInitiation(3),
      GetEnergy(kBulge1));
  EXPECT_EQ(em->InternalLoopInitiation(5) + std::min(em->internal_asym, NINIO_MAX_ASYM) +
          em->internal_au_penalty + em->au_penalty * 2 + em->internal_2x3_mismatch[A][G][A][U] +
          em->internal_2x3_mismatch[C][A][A][G] + em->HairpinInitiation(3),
      GetEnergy(kInternal1));
}

INSTANTIATE_TEST_SUITE_P(EnergyModelTests, T04ModelTest, testing::Range(0, NUM_TEST_MODELS));

struct CtdTest {
  Primary r;
  Secondary s;
  Ctds ctd;
  BranchCtd branch_ctd;
  std::deque<int> branches;
};

#if ENERGY_PRECISION == 1

TEST(T04P1ModelTest, T04P1) {
  auto em = t04p1;

  EXPECT_EQ(E(8.8), em->HairpinInitiation(87));
  EXPECT_EQ(E(6.8), em->BulgeInitiation(57));
  EXPECT_EQ(E(4.6), em->InternalLoopInitiation(67));

  EXPECT_EQ(E(4.5), GetEnergy(em, "GCAAAGCC", "((...).)"));
  EXPECT_EQ(E(5.7), GetEnergy(em, "CCCAAAAUG", ".(.(...))"));
  EXPECT_EQ(E(5.5), GetEnergy(em, "UACAGA", "(....)"));
  EXPECT_EQ(E(-0.6), GetEnergy(em, "AGGGUCAUCCG", ".(((...)))."));
  EXPECT_EQ(E(8.0), GetEnergy(em, "AGAGAAACAAAU", "(..(...)...)"));
  EXPECT_EQ(E(9.5), GetEnergy(em, "CGUUGCCUAAAAAGGAAACAAG", "(.............(...)..)"));
  EXPECT_EQ(E(7.7), GetEnergy(em, "CCCGAAACAG", "(..(...).)"));
  EXPECT_EQ(E(7.4), GetEnergy(em, "GACAGAAACGCUGAAUC", "((..(...)......))"));
  EXPECT_EQ(E(17.3), GetEnergy(em, "CUGAAACUGGAAACAGAAAUG", "(.(...)..(...).(...))"));
  EXPECT_EQ(E(18.2), GetEnergy(em, "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))"));
  EXPECT_EQ(E(17.6), GetEnergy(em, "AGCUAAAAACAAAGGUGAAACGU", "(..(...).(...)..(...).)"));
  EXPECT_EQ(E(13.1), GetEnergy(em, "CUGAAACUGGAAACAGAAAUG", ".(.(...)(....)......)"));
  EXPECT_EQ(E(-27.6),
      GetEnergy(em, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
          "(((((((((((.((...((((....))))..)).)))..((((..((((....))))...)))).))))))))...."));
  EXPECT_EQ(E(17.9), GetEnergy(em, "UCUGAGUAAAUUGCUACGCG", "(....)((...).......)"));
  EXPECT_EQ(E(-43.1), GetEnergy(em, k16sHSapiens3));

  // Special stacking - this is not implemented. TODO(4): Implement this?
  EXPECT_EQ(E(3.7), GetEnergy(em, "GGUCAAAGGUC", "((((...))))"));
  EXPECT_EQ(E(-4.5), GetEnergy(em, "GGGGAAACCCC", "((((...))))"));
  EXPECT_EQ(E(7.2), GetEnergy(em, "UGACAAAGGCGA", "(..(...)...)"));

  // NNDB flush coax
  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][A][U][G] +
          2 * em->au_penalty + 2 * em->HairpinInitiation(3),
      GetEnergy(em, "GUGAAACACAAAAUGA", ".((...))((...))."));

  // NNDB T99 Multiloop example
  EXPECT_EQ(em->stack[G][A][U][C] + em->terminal[C][G][A][G] + em->coax_mismatch_non_contiguous +
          3 * em->HairpinInitiation(3) + em->MultiloopInitiation(4) + 2 * em->au_penalty,
      GetEnergy(em, "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))"));
}

TEST(T04P1ModelTest, Precomp) {
  auto em = t04p1;

  const Precomp pc(Primary::FromSeq("GGGGAAACCCC"), em);
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

TEST(T04P2ModelTest, T04P2) {
  auto em = t04p2;

  EXPECT_EQ(E(8.85), em->HairpinInitiation(87));
  EXPECT_EQ(E(6.79), em->BulgeInitiation(57));
  EXPECT_EQ(E(4.57), em->InternalLoopInitiation(67));

  EXPECT_EQ(E(4.45), GetEnergy(em, "GCAAAGCC", "((...).)"));
  EXPECT_EQ(E(5.71), GetEnergy(em, "CCCAAAAUG", ".(.(...))"));
  EXPECT_EQ(E(5.50), GetEnergy(em, "UACAGA", "(....)"));
  EXPECT_EQ(E(-0.59), GetEnergy(em, "AGGGUCAUCCG", ".(((...)))."));
  EXPECT_EQ(E(8.00), GetEnergy(em, "AGAGAAACAAAU", "(..(...)...)"));
  EXPECT_EQ(E(9.50), GetEnergy(em, "CGUUGCCUAAAAAGGAAACAAG", "(.............(...)..)"));
  EXPECT_EQ(E(7.70), GetEnergy(em, "CCCGAAACAG", "(..(...).)"));
  EXPECT_EQ(E(7.45), GetEnergy(em, "GACAGAAACGCUGAAUC", "((..(...)......))"));
  EXPECT_EQ(E(17.30), GetEnergy(em, "CUGAAACUGGAAACAGAAAUG", "(.(...)..(...).(...))"));
  EXPECT_EQ(E(18.25), GetEnergy(em, "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))"));
  EXPECT_EQ(E(17.60), GetEnergy(em, "AGCUAAAAACAAAGGUGAAACGU", "(..(...).(...)..(...).)"));
  EXPECT_EQ(E(13.10), GetEnergy(em, "CUGAAACUGGAAACAGAAAUG", ".(.(...)(....)......)"));
  EXPECT_EQ(E(-27.41),
      GetEnergy(em, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
          "(((((((((((.((...((((....))))..)).)))..((((..((((....))))...)))).))))))))...."));
  EXPECT_EQ(E(17.90), GetEnergy(em, "UCUGAGUAAAUUGCUACGCG", "(....)((...).......)"));
  EXPECT_EQ(E(-42.75), GetEnergy(em, k16sHSapiens3));

  // Special stacking - this is not implemented. TODO(4): Implement this?
  EXPECT_EQ(E(3.63), GetEnergy(em, "GGUCAAAGGUC", "((((...))))"));
  EXPECT_EQ(E(-4.38), GetEnergy(em, "GGGGAAACCCC", "((((...))))"));
  EXPECT_EQ(E(7.20), GetEnergy(em, "UGACAAAGGCGA", "(..(...)...)"));

  // NNDB flush coax
  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][A][U][G] +
          2 * em->au_penalty + 2 * em->HairpinInitiation(3),
      GetEnergy(em, "GUGAAACACAAAAUGA", ".((...))((...))."));

  // NNDB T99 Multiloop example
  EXPECT_EQ(em->stack[G][A][U][C] + em->terminal[C][G][A][G] + em->coax_mismatch_non_contiguous +
          3 * em->HairpinInitiation(3) + em->MultiloopInitiation(4) + 2 * em->au_penalty,
      GetEnergy(em, "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))"));
}

TEST(T04P2ModelTest, Precomp) {
  auto em = t04p2;

  const Precomp pc(Primary::FromSeq("GGGGAAACCCC"), em);
  EXPECT_EQ(E(-2.10 - 0.40 - 1.60), pc.min_mismatch_coax);
  EXPECT_EQ(E(-3.42), pc.min_flush_coax);
  EXPECT_EQ(E(-2.60), pc.min_twoloop_not_stack);

  Energy augubranch[4][4] = {{E(-0.60), E(-0.60), E(-0.60), E(0.50 - 0.60)},
      {E(-0.60), E(-0.60), E(-0.60), E(-0.60)}, {E(-0.60), E(-0.60), E(-0.60), E(0.50 - 0.60)},
      {E(0.50 - 0.60), E(-0.60), E(0.50 - 0.60), E(-0.60)}};
  EXPECT_EQ(sizeof(augubranch), sizeof(pc.augubranch));
  EXPECT_EQ(0, std::memcmp(augubranch, pc.augubranch, sizeof(augubranch)));
}

TEST(T04P2ModelTest, T12P2) {
  auto em = t12p2;

  EXPECT_EQ(E(8.85), em->HairpinInitiation(87));
  EXPECT_EQ(E(6.79), em->BulgeInitiation(57));
  EXPECT_EQ(E(4.57), em->InternalLoopInitiation(67));

  // Example from https://doi.org/10.1093/nar/gkac261
  EXPECT_EQ(E(-1.74), GetEnergy(em, "UGUCGAUACCCUGUCGAUA", "((((((((...))))))))"));
  EXPECT_EQ(E(-3.08), GetEnergy(em, "UAGGUCAGCCCCUGGUCUA", "((((((((...))))))))"));

  EXPECT_EQ(E(4.45), GetEnergy(em, "GCAAAGCC", "((...).)"));
  EXPECT_EQ(E(5.71), GetEnergy(em, "CCCAAAAUG", ".(.(...))"));
  EXPECT_EQ(E(5.50), GetEnergy(em, "UACAGA", "(....)"));
  EXPECT_EQ(E(-1.36), GetEnergy(em, "AGGGUCAUCCG", ".(((...)))."));
  EXPECT_EQ(E(8.00), GetEnergy(em, "AGAGAAACAAAU", "(..(...)...)"));
  EXPECT_EQ(E(9.50), GetEnergy(em, "CGUUGCCUAAAAAGGAAACAAG", "(.............(...)..)"));
  EXPECT_EQ(E(7.70), GetEnergy(em, "CCCGAAACAG", "(..(...).)"));
  EXPECT_EQ(E(7.45), GetEnergy(em, "GACAGAAACGCUGAAUC", "((..(...)......))"));
  EXPECT_EQ(E(16.30), GetEnergy(em, "CUGAAACUGGAAACAGAAAUG", "(.(...)..(...).(...))"));
  EXPECT_EQ(E(18.25), GetEnergy(em, "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))"));
  EXPECT_EQ(E(17.60), GetEnergy(em, "AGCUAAAAACAAAGGUGAAACGU", "(..(...).(...)..(...).)"));
  EXPECT_EQ(E(12.10), GetEnergy(em, "CUGAAACUGGAAACAGAAAUG", ".(.(...)(....)......)"));
  EXPECT_EQ(E(-28.29),
      GetEnergy(em, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA",
          "(((((((((((.((...((((....))))..)).)))..((((..((((....))))...)))).))))))))...."));
  EXPECT_EQ(E(15.90), GetEnergy(em, "UCUGAGUAAAUUGCUACGCG", "(....)((...).......)"));
  EXPECT_EQ(E(-41.46), GetEnergy(em, k16sHSapiens3));

  // Special stacking - this is not implemented. TODO(4): Implement this?
  EXPECT_EQ(E(2.52), GetEnergy(em, "GGUCAAAGGUC", "((((...))))"));
  EXPECT_EQ(E(-4.38), GetEnergy(em, "GGGGAAACCCC", "((((...))))"));
  EXPECT_EQ(E(7.20), GetEnergy(em, "UGACAAAGGCGA", "(..(...)...)"));

  // NNDB flush coax
  EXPECT_EQ(em->stack[C][A][U][G] + em->stack[A][C][G][U] + em->stack[C][A][U][G] +
          2 * em->au_penalty + 2 * em->HairpinInitiation(3),
      GetEnergy(em, "GUGAAACACAAAAUGA", ".((...))((...))."));

  // NNDB T99 Multiloop example
  EXPECT_EQ(em->stack[G][A][U][C] + em->terminal[C][G][A][G] + em->coax_mismatch_non_contiguous +
          3 * em->HairpinInitiation(3) + em->MultiloopInitiation(4) + 2 * em->au_penalty,
      GetEnergy(em, "UUAGAAACGCAAAGAGGUCCAAAGA", "(..(...).(...).....(...))"));
}

// NEWMODEL: Add tests here.

#endif

std::function<CtdTest(const Model::Ptr&)> CTD_TESTS[] = {[](const Model::Ptr&) -> CtdTest {
                                                           return {{}, {}, {}, {}, {}};
                                                         },
    [](const Model::Ptr&) -> CtdTest {
      return {Primary::FromSeq("A"), Secondary::FromDb("."), Ctds{CTD_NA}, {}, {}};
    },
    [](const Model::Ptr&) -> CtdTest {
      return {Primary::FromSeq("AG"), Secondary::FromDb(".."), Ctds{CTD_NA, CTD_NA}, {}, {}};
    },
    [](const Model::Ptr&) -> CtdTest {
      return {
          Primary::FromSeq("GUA"), Secondary::FromDb("..."), Ctds{CTD_NA, CTD_NA, CTD_NA}, {}, {}};
    },
    [](const Model::Ptr&) -> CtdTest {
      return {Primary::FromSeq("GUAC"), Secondary::FromDb("...."),
          Ctds{CTD_NA, CTD_NA, CTD_NA, CTD_NA}, {}, {}};
    },
    // 3' dangle inside the branch.
    [](const Model::Ptr& em) -> CtdTest {
      return {Primary::FromSeq("GAAAC"), Secondary::FromDb("(...)"),
          Ctds{CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_3_DANGLE},
          {{CTD_3_DANGLE, em->dangle3[G][A][C]}}, {4}};
    },
    [](const Model::Ptr&) -> CtdTest {
      return {Primary::FromSeq("GAAACAGAAAAUGGAAACCAGAAACA"),
          Secondary::FromDb("(...).((...).(...)).(...)."), Ctds(26), {}, {}};
    },
    [](const Model::Ptr& em) -> CtdTest {
      return {Primary::FromSeq("GAAACAGAAAAUGGAAACCAGAAACA"),
          Secondary::FromDb("(...).((...).(...)).(...)."),
          Ctds{CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_RC_WITH_NEXT, CTD_NA, CTD_NA,
              CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
              CTD_NA, CTD_RC_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA},
          {{CTD_UNUSED, ZERO_E}, {CTD_RC_WITH_NEXT, em->MismatchCoaxial(C, A, A, G)},
              {CTD_RC_WITH_PREV, em->MismatchCoaxial(C, A, A, G)}},
          {0, 6, 20}};
    },
    [](const Model::Ptr& em) -> CtdTest {
      return {Primary::FromSeq("GAAACAGAAAAUGGAAACCAGAAACA"),
          Secondary::FromDb("(...).((...).(...)).(...)."),
          Ctds{CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_PREV, CTD_NA,
              CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_5_DANGLE, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
              CTD_FCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA},
          {{CTD_FCOAX_WITH_NEXT, em->stack[G][A][U][C]},
              {CTD_FCOAX_WITH_PREV, em->stack[G][A][U][C]}, {CTD_5_DANGLE, em->dangle5[C][G][G]}},
          {18, 7, 13}};
    },
    [](const Model::Ptr& em) -> CtdTest {
      return {Primary::FromSeq("GGAAACGAAACC"), Secondary::FromDb("((...)(...))"),
          Ctds{CTD_NA, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_NEXT, CTD_NA,
              CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_PREV},
          {{CTD_UNUSED, ZERO_E}, {CTD_FCOAX_WITH_NEXT, em->stack[G][G][C][C]},
              {CTD_FCOAX_WITH_PREV, em->stack[G][G][C][C]}},
          {1, 6, 11}};
    },
    [](const Model::Ptr& em) -> CtdTest {
      return {Primary::FromSeq("UUAGAAACGCAAAGAGGUCCAAAGA"),
          Secondary::FromDb("(..(...).(...).....(...))"),
          Ctds{CTD_NA, CTD_NA, CTD_NA, CTD_LCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
              CTD_LCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
              CTD_NA, CTD_FCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_PREV},
          {{CTD_FCOAX_WITH_PREV, em->stack[U][C][G][A]},
              {CTD_LCOAX_WITH_NEXT, em->MismatchCoaxial(C, G, A, G)},
              {CTD_LCOAX_WITH_PREV, em->MismatchCoaxial(C, G, A, G)},
              {CTD_FCOAX_WITH_NEXT, em->stack[U][C][G][A]}},
          {24, 3, 9, 19}};
    }};

class CtdsTest
    : public testing::TestWithParam<std::tuple<int, std::function<CtdTest(const Model::Ptr&)>>> {};

TEST_P(CtdsTest, BaseBranchBase) {
  auto em = test_t04_ems[std::get<0>(GetParam())];
  auto ctd_test = std::get<1>(GetParam())(em);
  // Convert base representation to branch representation.
  BranchCtd computed_branch_ctd;
  auto computed_energy = AddBaseCtdsToBranchCtds(
      *em, ctd_test.r, ctd_test.s, ctd_test.ctd, ctd_test.branches, &computed_branch_ctd);
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

INSTANTIATE_TEST_SUITE_P(CtdsTest, CtdsTest,
    testing::Combine(testing::Range(0, NUM_TEST_MODELS), testing::ValuesIn(CTD_TESTS)));

}  // namespace mrna::md::t04
