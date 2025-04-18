// Copyright 2024 E.
#include "api/pfn.h"

#include <string>
#include <tuple>

#include "api/ctx/ctx_cfg.h"
#include "gtest/gtest.h"
#include "model/primary.h"
#include "tests/init.h"
#include "tests/util.h"
#include "util/float.h"

namespace mrna {

class PfnTestT04 : public testing::TestWithParam<std::tuple<int, CtxCfg::PfnAlg>> {
 public:
  static pfn::PfnResult Pfn(const BackendModelPtr& m, const std::string& s) {
    return GetPfn(m, std::get<1>(GetParam()), s);
  }

  static pfn::PfnResult Pfn(const BackendModelPtr& m, const Primary& r) {
    return GetPfn(m, std::get<1>(GetParam()), r);
  }
};

#if ENERGY_PRECISION == 1

TEST_P(PfnTestT04, T04P1) {
  auto [i, alg] = GetParam();
  auto m = t04_ms[i];
  if (!Contains(CtxCfg::PfnAlgsForBackend(m), alg)) return;

  EXPECT_REL_EQ(FLT(4.2481601382949495665565296828679689667765375832), Pfn(m, "CCUCCGGG").pfn.q);
  EXPECT_REL_EQ(FLT(4.17979557041608366287852107192666645517291810433), Pfn(m, "CGGAAACGG").pfn.q);
  EXPECT_REL_EQ(FLT(4.62750397465096928509866459188973500341174374984), Pfn(m, "UGCAAAGCAA").pfn.q);
  EXPECT_REL_EQ(
      FLT(3122.66115008916122674423382261602803687854854357), Pfn(m, "GGGGAAACCCC").pfn.q);
  EXPECT_REL_EQ(
      FLT(9.44718000805199118395269282804335660230798127861), Pfn(m, "CUUAUAGUUAAGG").pfn.q);
  EXPECT_REL_EQ(
      FLT(930.051311983052209666709379864276379348510979399), Pfn(m, "CCGAAGGGGCUGCGGCG").pfn.q);
  EXPECT_REL_EQ(
      FLT(196.739899801268074338798758309195750804501201612), Pfn(m, "GCCAAGGCCCCACCCGGA").pfn.q);
  EXPECT_REL_EQ(
      FLT(3431.19128512585938334129785636062228774823050418), Pfn(m, "GGCCGAUGGCAGCGAUAGC").pfn.q);
  EXPECT_REL_EQ(FLT(94.3008892348264464511025477744977107828479729955),
      Pfn(m, "CUGAAACUGGAAACAGAAAUG").pfn.q);

  // Too slow for brute force:
  if (alg == CtxCfg::PfnAlg::BRUTE) return;
  EXPECT_REL_EQ(FLT(573963557.832690101314804611565712243969839876608),
      Pfn(m, "CCGGGCCAGCCCGCUCCUACGGGGGGUC").pfn.q);
  EXPECT_REL_EQ(FLT(226979.219096921172665125564673329771324145408152),
      Pfn(m, "CGCAGGGUCGGACCCGGGAGAACCGCGA").pfn.q);
  EXPECT_REL_EQ(FLT(812.383523749713349970768190632670001171883879147),
      Pfn(m, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG").pfn.q);
  EXPECT_REL_EQ(FLT(19520.3202391172281998883152194957928196468529852),
      Pfn(m, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU").pfn.q);
  EXPECT_REL_EQ(FLT(2003.03838946282814030479317865110152744930660671),
      Pfn(m, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG").pfn.q);
  EXPECT_REL_EQ(FLT(23058.3336091435782209696760380319112512888170753),
      Pfn(m, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC").pfn.q);
  EXPECT_REL_EQ(FLT(2060.71972602168776467337879240845814018091207899),
      Pfn(m, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC").pfn.q);
  EXPECT_REL_EQ(FLT(10731115.7023640924024241047184987959200656768803),
      Pfn(m, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG").pfn.q);
  EXPECT_REL_EQ(FLT(443426.84236219495995929904990764754412943693663),
      Pfn(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC").pfn.q);
  EXPECT_REL_EQ(FLT(349171.530355094625444533151609745244789447612351),
      Pfn(m, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA").pfn.q);
  EXPECT_REL_EQ(FLT(132689.019044556375490017850862163999027602957435),
      Pfn(m, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC").pfn.q);
  EXPECT_REL_EQ(FLT(803546.641992943133871048721772712844519695820412),
      Pfn(m, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG").pfn.q);
  EXPECT_REL_EQ(FLT(2520.43161740232401829233973384101843671743563485),
      Pfn(m, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG").pfn.q);
  EXPECT_REL_EQ(FLT(643578.928228185905419648713152576603855800599177),
      Pfn(m, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG").pfn.q);
  EXPECT_REL_EQ(FLT(1660031715.98196937466554739068134244758781101395),
      Pfn(m, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG").pfn.q);
  EXPECT_REL_EQ(FLT(17711.6389272494574341953125850175769175966289746),
      Pfn(m, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC").pfn.q);
  EXPECT_REL_EQ(FLT(413455851608.007014854697696998395303922733033517),
      Pfn(m, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG").pfn.q);
  EXPECT_REL_EQ(FLT(25735350509881543.1432642249111720116326051872321),
      Pfn(m, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU").pfn.q);
  EXPECT_REL_EQ(FLT(265269.030237070990473213695638615156718492406518),
      Pfn(m, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA").pfn.q);
  EXPECT_REL_EQ(FLT(19869500930.8524243679823577324638360720604214023),
      Pfn(m, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC").pfn.q);
  EXPECT_REL_EQ(FLT(355492874401.288593060501840506007538647116228473),
      Pfn(m, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU").pfn.q);
  EXPECT_REL_EQ(FLT(1404779993264957699657.30024473387442890744595273),
      Pfn(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA")
          .pfn.q);
  EXPECT_REL_EQ(FLT(7.19185047044768251505176115957271369407798113857e+68),
      Pfn(m, std::get<Primary>(k16sHSapiens3)).pfn.q);
}

#elif ENERGY_PRECISION == 2

TEST_P(PfnTestT04, T04P2) {
  auto [i, alg] = GetParam();
  const auto& m = t04_ms[i];
  if (!Contains(CtxCfg::PfnAlgsForBackend(m), alg)) return;

  // Regression tests:
  auto pfn = Pfn(m, "GGCGACCGGCGG").pfn;
  PfnTables res(BoltzSums(12, FLT(0.0)), FLT(40.3707273271976), BoltzProbs(12, FLT(0.0)));
  res.p[0][6] = FLT(0.03214767742222614);
  res.p[1][5] = FLT(0.0001566462984081384);
  res.p[1][9] = FLT(1.215299560015079);
  res.p[2][8] = FLT(0.004727955899485569);
  res.p[2][10] = FLT(0.03719423081290614);
  res.p[3][9] = FLT(0.0007935322030592623);
  res.p[5][1] = FLT(3335.3194400046264);
  res.p[5][11] = FLT(0.03153767346064938);
  res.p[6][0] = FLT(16.772317556168453);
  res.p[6][10] = FLT(0.0001566462984081384);
  res.p[8][2] = FLT(8023.220887240605);
  res.p[9][1] = FLT(31.182388349931827);
  res.p[9][3] = FLT(964.9312203906061);
  res.p[10][2] = FLT(20.319418786068947);
  res.p[10][6] = FLT(666.6532747397315);
  res.p[11][5] = FLT(3.2507238321139087);
  res.prob[0][6] = FLT(0.013355990593104645);
  res.prob[1][5] = FLT(0.012941690152147614);
  res.prob[1][9] = FLT(0.9386985410183861);
  res.prob[2][8] = FLT(0.9396272259169692);
  res.prob[2][10] = FLT(0.018720622647885602);
  res.prob[3][9] = FLT(0.018966812039607928);
  res.prob[5][11] = FLT(0.002539470391431164);
  res.prob[6][10] = FLT(0.002586744770864938);

  CheckPfn(pfn, res);

  EXPECT_REL_EQ(FLT(4.0626572470100567109913456873125596215663636765), Pfn(m, "CCUCCGGG").pfn.q);
  EXPECT_REL_EQ(FLT(3.99326566300301791033216574191242303938854885947), Pfn(m, "CGGAAACGG").pfn.q);
  EXPECT_REL_EQ(FLT(4.78710233484875966364273595867119395656372672638), Pfn(m, "UGCAAAGCAA").pfn.q);
  EXPECT_REL_EQ(
      FLT(2662.62888037900558704265668051982251217213153673), Pfn(m, "GGGGAAACCCC").pfn.q);
  EXPECT_REL_EQ(
      FLT(9.99242559078143218348552048735575294737368979132), Pfn(m, "CUUAUAGUUAAGG").pfn.q);
  EXPECT_REL_EQ(
      FLT(862.833647006987128073360487995184990796940938988), Pfn(m, "CCGAAGGGGCUGCGGCG").pfn.q);
  EXPECT_REL_EQ(
      FLT(187.240327086207995339175981065221999797886386151), Pfn(m, "GCCAAGGCCCCACCCGGA").pfn.q);
  EXPECT_REL_EQ(
      FLT(3426.95720175447283008095614881832914527009516235), Pfn(m, "GGCCGAUGGCAGCGAUAGC").pfn.q);
  EXPECT_REL_EQ(FLT(92.8392679166754970004496359132818709013811507588),
      Pfn(m, "CUGAAACUGGAAACAGAAAUG").pfn.q);

  // Too slow for brute force:
  if (alg == CtxCfg::PfnAlg::BRUTE) return;
  EXPECT_REL_EQ(FLT(485598279.6070418149352955163452483973488643105),
      Pfn(m, "CCGGGCCAGCCCGCUCCUACGGGGGGUC").pfn.q);
  EXPECT_REL_EQ(FLT(185629.387202844402235328214859370318079740293414),
      Pfn(m, "CGCAGGGUCGGACCCGGGAGAACCGCGA").pfn.q);
  EXPECT_REL_EQ(FLT(667.550167100069406093819365050392678255967576263),
      Pfn(m, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG").pfn.q);
  EXPECT_REL_EQ(FLT(19052.940428268950747234584434128107007484910591),
      Pfn(m, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU").pfn.q);
  EXPECT_REL_EQ(FLT(1801.93385140750909623037431133829103418741017157),
      Pfn(m, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG").pfn.q);
  EXPECT_REL_EQ(FLT(23323.1953589110134469979245648405509649338667612),
      Pfn(m, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC").pfn.q);
  EXPECT_REL_EQ(FLT(1851.60443772924479733223300902619145544021130277),
      Pfn(m, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC").pfn.q);
  EXPECT_REL_EQ(FLT(11249868.4390407380833978400336731868612605274107),
      Pfn(m, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG").pfn.q);
  EXPECT_REL_EQ(FLT(435607.642597553412901094152749505662052383450608),
      Pfn(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC").pfn.q);
  EXPECT_REL_EQ(FLT(322009.249593371879816329736817237306336628900648),
      Pfn(m, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA").pfn.q);
  EXPECT_REL_EQ(FLT(118810.21028108107308050449782865820084230958938),
      Pfn(m, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC").pfn.q);
  EXPECT_REL_EQ(FLT(749592.567993534559924070154970812994457327672832),
      Pfn(m, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG").pfn.q);
  EXPECT_REL_EQ(FLT(2563.33480855576407734846090041277273681152434772),
      Pfn(m, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG").pfn.q);
  EXPECT_REL_EQ(FLT(565652.958252855139999190467774364353485165171175),
      Pfn(m, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG").pfn.q);
  EXPECT_REL_EQ(FLT(1481396422.98096632115278359427469196950840412232),
      Pfn(m, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG").pfn.q);
  EXPECT_REL_EQ(FLT(18301.1241716188359533542291829161488228425372024),
      Pfn(m, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC").pfn.q);
  EXPECT_REL_EQ(FLT(384503084677.407285224649566220747528595879119908),
      Pfn(m, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG").pfn.q);
  EXPECT_REL_EQ(FLT(18644551213478601.4954775128877259814283674980883),
      Pfn(m, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU").pfn.q);
  EXPECT_REL_EQ(FLT(237913.207098452661898475512667249518513054819729),
      Pfn(m, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA").pfn.q);
  EXPECT_REL_EQ(FLT(15152854806.8121027817982722467273280176326700692),
      Pfn(m, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC").pfn.q);
  EXPECT_REL_EQ(FLT(311680103014.139072868971361273691609978660113255),
      Pfn(m, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU").pfn.q);
  EXPECT_REL_EQ(FLT(855252103666168038753.963564419273686450258512786),
      Pfn(m, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA")
          .pfn.q);
  EXPECT_REL_EQ(FLT(4.03150179479028914802269005329388617139041407731e+68),
      Pfn(m, std::get<Primary>(k16sHSapiens3)).pfn.q);
}

#endif

INSTANTIATE_TEST_SUITE_P(PfnTest, PfnTestT04,
    testing::Combine(
        testing::Range(0, NUM_T04_MODELS), testing::ValuesIn(EnumValues<CtxCfg::PfnAlg>())));

}  // namespace mrna
