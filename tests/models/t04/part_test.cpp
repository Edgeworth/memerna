// Copyright 2022 E.
#include "api/part.h"

#include <string>
#include <tuple>

#include "api/ctx/ctx.h"
#include "api/ctx/ctx_cfg.h"
#include "api/energy/model.h"
#include "common_test.h"
#include "gtest/gtest.h"
#include "model/primary.h"
#include "util/float.h"

namespace mrna::md::t04 {

class PartitionTestT04Like : public testing::TestWithParam<CtxCfg::PartAlg> {
 public:
  static part::PartResult Partition(const erg::EnergyModelPtr& em, const std::string& s) {
    return Partition(em, Primary::FromSeq(s));
  }

  static part::PartResult Partition(const erg::EnergyModelPtr& em, const Primary& r) {
    return Ctx(em, CtxCfg{.part_alg = GetParam()}).Partition(r);
  }
};

#if ENERGY_PRECISION == 1

TEST_P(PartitionTestT04Like, T04P1) {
  auto em = t04p1;

  EXPECT_REL_EQ(
      FLT(4.2481601382949495665565296828679689667765375832), Partition(em, "CCUCCGGG").part.q);
  EXPECT_REL_EQ(
      FLT(4.17979557041608366287852107192666645517291810433), Partition(em, "CGGAAACGG").part.q);
  EXPECT_REL_EQ(
      FLT(4.62750397465096928509866459188973500341174374984), Partition(em, "UGCAAAGCAA").part.q);
  EXPECT_REL_EQ(
      FLT(3122.66115008916122674423382261602803687854854357), Partition(em, "GGGGAAACCCC").part.q);
  EXPECT_REL_EQ(FLT(9.44718000805199118395269282804335660230798127861),
      Partition(em, "CUUAUAGUUAAGG").part.q);
  EXPECT_REL_EQ(FLT(930.051311983052209666709379864276379348510979399),
      Partition(em, "CCGAAGGGGCUGCGGCG").part.q);
  EXPECT_REL_EQ(FLT(196.739899801268074338798758309195750804501201612),
      Partition(em, "GCCAAGGCCCCACCCGGA").part.q);
  EXPECT_REL_EQ(FLT(3431.19128512585938334129785636062228774823050418),
      Partition(em, "GGCCGAUGGCAGCGAUAGC").part.q);
  EXPECT_REL_EQ(FLT(94.3008892348264464511025477744977107828479729955),
      Partition(em, "CUGAAACUGGAAACAGAAAUG").part.q);

  // Too slow for brute force:
  if (GetParam() == CtxCfg::PartAlg::BRUTE) return;
  EXPECT_REL_EQ(FLT(573963557.832690101314804611565712243969839876608),
      Partition(em, "CCGGGCCAGCCCGCUCCUACGGGGGGUC").part.q);
  EXPECT_REL_EQ(FLT(226979.219096921172665125564673329771324145408152),
      Partition(em, "CGCAGGGUCGGACCCGGGAGAACCGCGA").part.q);
  EXPECT_REL_EQ(FLT(812.383523749713349970768190632670001171883879147),
      Partition(em, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG").part.q);
  EXPECT_REL_EQ(FLT(19520.3202391172281998883152194957928196468529852),
      Partition(em, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU").part.q);
  EXPECT_REL_EQ(FLT(2003.03838946282814030479317865110152744930660671),
      Partition(em, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG").part.q);
  EXPECT_REL_EQ(FLT(23058.3336091435782209696760380319112512888170753),
      Partition(em, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC").part.q);
  EXPECT_REL_EQ(FLT(2060.71972602168776467337879240845814018091207899),
      Partition(em, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC").part.q);
  EXPECT_REL_EQ(FLT(10731115.7023640924024241047184987959200656768803),
      Partition(em, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG").part.q);
  EXPECT_REL_EQ(FLT(443426.84236219495995929904990764754412943693663),
      Partition(em, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC").part.q);
  EXPECT_REL_EQ(FLT(349171.530355094625444533151609745244789447612351),
      Partition(em, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA").part.q);
  EXPECT_REL_EQ(FLT(132689.019044556375490017850862163999027602957435),
      Partition(em, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC").part.q);
  EXPECT_REL_EQ(FLT(803546.641992943133871048721772712844519695820412),
      Partition(em, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG").part.q);
  EXPECT_REL_EQ(FLT(2520.43161740232401829233973384101843671743563485),
      Partition(em, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG").part.q);
  EXPECT_REL_EQ(FLT(643578.928228185905419648713152576603855800599177),
      Partition(em, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG").part.q);
  EXPECT_REL_EQ(FLT(1660031715.98196937466554739068134244758781101395),
      Partition(em, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG").part.q);
  EXPECT_REL_EQ(FLT(17711.6389272494574341953125850175769175966289746),
      Partition(em, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC").part.q);
  EXPECT_REL_EQ(FLT(413455851608.007014854697696998395303922733033517),
      Partition(em, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG").part.q);
  EXPECT_REL_EQ(FLT(25735350509881543.1432642249111720116326051872321),
      Partition(em, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU").part.q);
  EXPECT_REL_EQ(FLT(265269.030237070990473213695638615156718492406518),
      Partition(em, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA").part.q);
  EXPECT_REL_EQ(FLT(19869500930.8524243679823577324638360720604214023),
      Partition(em, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC").part.q);
  EXPECT_REL_EQ(FLT(355492874401.288593060501840506007538647116228473),
      Partition(em, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU")
          .part.q);
  EXPECT_REL_EQ(FLT(1404779993264957699657.30024473387442890744595273),
      Partition(em, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA")
          .part.q);
  EXPECT_REL_EQ(FLT(7.19185047044768251505176115957271369407798113857e+68),
      Partition(em, std::get<Primary>(k16sHSapiens3)).part.q);
}

#elif ENERGY_PRECISION == 2

TEST_P(PartitionTestT04Like, T04P2) {
  auto em = t04p2;

  EXPECT_REL_EQ(
      FLT(4.0626572470100567109913456873125596215663636765), Partition(em, "CCUCCGGG").part.q);
  EXPECT_REL_EQ(
      FLT(3.99326566300301791033216574191242303938854885947), Partition(em, "CGGAAACGG").part.q);
  EXPECT_REL_EQ(
      FLT(4.78710233484875966364273595867119395656372672638), Partition(em, "UGCAAAGCAA").part.q);
  EXPECT_REL_EQ(
      FLT(2662.62888037900558704265668051982251217213153673), Partition(em, "GGGGAAACCCC").part.q);
  EXPECT_REL_EQ(FLT(9.99242559078143218348552048735575294737368979132),
      Partition(em, "CUUAUAGUUAAGG").part.q);
  EXPECT_REL_EQ(FLT(862.833647006987128073360487995184990796940938988),
      Partition(em, "CCGAAGGGGCUGCGGCG").part.q);
  EXPECT_REL_EQ(FLT(187.240327086207995339175981065221999797886386151),
      Partition(em, "GCCAAGGCCCCACCCGGA").part.q);
  EXPECT_REL_EQ(FLT(3426.95720175447283008095614881832914527009516235),
      Partition(em, "GGCCGAUGGCAGCGAUAGC").part.q);
  EXPECT_REL_EQ(FLT(92.8392679166754970004496359132818709013811507588),
      Partition(em, "CUGAAACUGGAAACAGAAAUG").part.q);

  // Too slow for brute force:
  if (GetParam() == CtxCfg::PartAlg::BRUTE) return;
  EXPECT_REL_EQ(FLT(485598279.6070418149352955163452483973488643105),
      Partition(em, "CCGGGCCAGCCCGCUCCUACGGGGGGUC").part.q);
  EXPECT_REL_EQ(FLT(185629.387202844402235328214859370318079740293414),
      Partition(em, "CGCAGGGUCGGACCCGGGAGAACCGCGA").part.q);
  EXPECT_REL_EQ(FLT(667.550167100069406093819365050392678255967576263),
      Partition(em, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG").part.q);
  EXPECT_REL_EQ(FLT(19052.940428268950747234584434128107007484910591),
      Partition(em, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU").part.q);
  EXPECT_REL_EQ(FLT(1801.93385140750909623037431133829103418741017157),
      Partition(em, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG").part.q);
  EXPECT_REL_EQ(FLT(23323.1953589110134469979245648405509649338667612),
      Partition(em, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC").part.q);
  EXPECT_REL_EQ(FLT(1851.60443772924479733223300902619145544021130277),
      Partition(em, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC").part.q);
  EXPECT_REL_EQ(FLT(11249868.4390407380833978400336731868612605274107),
      Partition(em, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG").part.q);
  EXPECT_REL_EQ(FLT(435607.642597553412901094152749505662052383450608),
      Partition(em, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC").part.q);
  EXPECT_REL_EQ(FLT(322009.249593371879816329736817237306336628900648),
      Partition(em, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA").part.q);
  EXPECT_REL_EQ(FLT(118810.21028108107308050449782865820084230958938),
      Partition(em, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC").part.q);
  EXPECT_REL_EQ(FLT(749592.567993534559924070154970812994457327672832),
      Partition(em, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG").part.q);
  EXPECT_REL_EQ(FLT(2563.33480855576407734846090041277273681152434772),
      Partition(em, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG").part.q);
  EXPECT_REL_EQ(FLT(565652.958252855139999190467774364353485165171175),
      Partition(em, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG").part.q);
  EXPECT_REL_EQ(FLT(1481396422.98096632115278359427469196950840412232),
      Partition(em, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG").part.q);
  EXPECT_REL_EQ(FLT(18301.1241716188359533542291829161488228425372024),
      Partition(em, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC").part.q);
  EXPECT_REL_EQ(FLT(384503084677.407285224649566220747528595879119908),
      Partition(em, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG").part.q);
  EXPECT_REL_EQ(FLT(18644551213478601.4954775128877259814283674980883),
      Partition(em, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU").part.q);
  EXPECT_REL_EQ(FLT(237913.207098452661898475512667249518513054819729),
      Partition(em, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA").part.q);
  EXPECT_REL_EQ(FLT(15152854806.8121027817982722467273280176326700692),
      Partition(em, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC").part.q);
  EXPECT_REL_EQ(FLT(311680103014.139072868971361273691609978660113255),
      Partition(em, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU")
          .part.q);
  EXPECT_REL_EQ(FLT(855252103666168038753.963564419273686450258512786),
      Partition(em, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA")
          .part.q);
  EXPECT_REL_EQ(FLT(4.03150179479028914802269005329388617139041407731e+68),
      Partition(em, std::get<Primary>(k16sHSapiens3)).part.q);
}

TEST_P(PartitionTestT04Like, T12P2) {
  auto em = t12p2;

  EXPECT_REL_EQ(
      FLT(4.06569633368939134129842956017929466372263904011), Partition(em, "CCUCCGGG").part.q);
  EXPECT_REL_EQ(
      FLT(3.99326566300301791033216574191242303938854885947), Partition(em, "CGGAAACGG").part.q);
  EXPECT_REL_EQ(
      FLT(4.78710233484875966364273595867119395656372672638), Partition(em, "UGCAAAGCAA").part.q);
  EXPECT_REL_EQ(
      FLT(2662.62888037900558704265668051982251217213153673), Partition(em, "GGGGAAACCCC").part.q);
  EXPECT_REL_EQ(FLT(9.99563943739342085541269179482384129892368533722),
      Partition(em, "CUUAUAGUUAAGG").part.q);
  EXPECT_REL_EQ(FLT(862.868063685304486819230022354325051191500088491),
      Partition(em, "CCGAAGGGGCUGCGGCG").part.q);
  EXPECT_REL_EQ(FLT(187.240327086207995339175981065221999797886386151),
      Partition(em, "GCCAAGGCCCCACCCGGA").part.q);
  EXPECT_REL_EQ(FLT(3427.25182440763216480351674919522989458198248782),
      Partition(em, "GGCCGAUGGCAGCGAUAGC").part.q);
  EXPECT_REL_EQ(FLT(95.3203576551913841120607445389624915106981460066),
      Partition(em, "CUGAAACUGGAAACAGAAAUG").part.q);

  // Too slow for brute force:
  if (GetParam() == CtxCfg::PartAlg::BRUTE) return;
  EXPECT_REL_EQ(FLT(489568270.164770419074506176947118047561497788815),
      Partition(em, "CCGGGCCAGCCCGCUCCUACGGGGGGUC").part.q);
  EXPECT_REL_EQ(FLT(187180.833955420153362511796970009127830070133408),
      Partition(em, "CGCAGGGUCGGACCCGGGAGAACCGCGA").part.q);
  EXPECT_REL_EQ(FLT(814.008918620087370568179693189208017088096864986),
      Partition(em, "UACCCUGUUCAGCAUUGGAAAUUUCCUGGG").part.q);
  EXPECT_REL_EQ(FLT(22741.5868387155689286279036191520339838073439982),
      Partition(em, "GCGCCCCAGUCGACGCUGAGCUCCUCUGCU").part.q);
  EXPECT_REL_EQ(FLT(1553.20048228785844588065846731295769815746227141),
      Partition(em, "CCCAACGGAGUAACUUAGCGAAUAGCAGGGG").part.q);
  EXPECT_REL_EQ(FLT(16506.5628374857920600909881467801681421370345844),
      Partition(em, "GGCGCACGCGUUAGCCGGGGAUCCACAGUGC").part.q);
  EXPECT_REL_EQ(FLT(2766.18219170267040968134969429655986933697602129),
      Partition(em, "CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC").part.q);
  EXPECT_REL_EQ(FLT(16314470.3366908442969861810014982898798240233482),
      Partition(em, "UCCACGGCUCGACGGCGCACUUAGUGCGUGGG").part.q);
  EXPECT_REL_EQ(FLT(174354.136340798747237142734124689962155219671546),
      Partition(em, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC").part.q);
  EXPECT_REL_EQ(FLT(306757.35359762952567953542324994395947672765144),
      Partition(em, "CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA").part.q);
  EXPECT_REL_EQ(FLT(124440.050413952833453556719729744722178663147699),
      Partition(em, "AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC").part.q);
  EXPECT_REL_EQ(FLT(927455.616926851270597371460103197308424506866451),
      Partition(em, "CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG").part.q);
  EXPECT_REL_EQ(FLT(3587.34358858919765684031015543333796302415958089),
      Partition(em, "GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG").part.q);
  EXPECT_REL_EQ(FLT(600608.359736951624049941754677153894661712659679),
      Partition(em, "CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG").part.q);
  EXPECT_REL_EQ(FLT(1481658714.9941444315310894152364794441233085645),
      Partition(em, "ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG").part.q);
  EXPECT_REL_EQ(FLT(15550.8855514817640395248766652356206702651430645),
      Partition(em, "GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC").part.q);
  EXPECT_REL_EQ(FLT(176566679168.85927333626088323239116642854614832),
      Partition(em, "GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG").part.q);
  EXPECT_REL_EQ(FLT(4659356908563728.47058108720744252481389275608726),
      Partition(em, "UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU").part.q);
  EXPECT_REL_EQ(FLT(307224.227055583795926137214729598592857602827188),
      Partition(em, "UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA").part.q);
  EXPECT_REL_EQ(FLT(24442100708.1476972328721860799593888394102495689),
      Partition(em, "GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC").part.q);
  EXPECT_REL_EQ(FLT(108845941078.995157523709039360707789300153288793),
      Partition(em, "AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU")
          .part.q);
  EXPECT_REL_EQ(FLT(5288640267480376963639.59972345718964029669389315),
      Partition(em, "GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA")
          .part.q);
  EXPECT_REL_EQ(FLT(2.2656269188107328870825780514941627464024016996e+68),
      Partition(em, std::get<Primary>(k16sHSapiens3)).part.q);
}

// NEWMODEL: Add tests here.

#endif

INSTANTIATE_TEST_SUITE_P(PartitionTest, PartitionTestT04Like,
    testing::ValuesIn(CtxCfg::PartitionAlgsForModelKind(erg::ModelKind::T04_LIKE)));

}  // namespace mrna::md::t04
