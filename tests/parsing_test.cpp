// Copyright 2016 Eliot Courtney.
#include "gtest/gtest.h"
#include "model/ctd.h"

namespace mrna {

class ParsingTest : public testing::Test {
 public:
  const Primary kPrimary1 = {};
  const Secondary kSecondary1 = {};
  const Ctds kCtd1 = {};
  const std::string kCtdString1 = "";

  const Primary kPrimary2 = Primary::FromString("A");
  const Secondary kSecondary2 = {-1};
  const Ctds kCtd2 = {CTD_NA};
  const std::string kCtdString2 = ".";

  const Primary kPrimary3 = Primary::FromString("UACGUUGGUGCUUA");
  const Secondary kSecondary3 = {13, 5, -1, -1, -1, 1, -1, 11, -1, -1, -1, 7, -1, 0};
  const Ctds kCtd3 = {CTD_UNUSED, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_3_DANGLE,
      CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_UNUSED};
  const std::string kCtdString3 = "[[...].[...]3]";

  const Primary kPrimary4 = Primary::FromString("UACGUUGGUGCUUA");
  const Secondary kSecondary4 = {-1, 5, -1, -1, -1, 1, -1, 11, -1, -1, -1, 7, -1, -1};
  const Ctds kCtd4 = {CTD_NA, CTD_RCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
      CTD_RCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA};
  const std::string kCtdString4 = ".n...]mp...]M.";

  const Primary kPrimary5 = Primary::FromString("UACGUUGGUGCUUA");
  const Secondary kSecondary5 = {-1, 5, -1, -1, -1, 1, -1, 11, -1, -1, -1, 7, -1, -1};
  const Ctds kCtd5 = {CTD_NA, CTD_3_DANGLE, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_UNUSED,
      CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA};
  const std::string kCtdString5 = ".[...]3[...]..";

  const Primary kPrimary6 = Primary::FromString("UACGUUGGUGCUUA");
  const Secondary kSecondary6 = {6, 5, -1, -1, -1, 1, 0, 11, -1, -1, -1, 7, -1, -1};
  const Ctds kCtd6 = {CTD_FCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
      CTD_FCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA};
  const std::string kCtdString6 = "n[...]]p...]..";

  const Primary kPrimary7 = Primary::FromString("UACGUUGGUGCUU");
  const Secondary kSecondary7 = {12, 5, -1, -1, -1, 1, -1, 11, -1, -1, -1, 7, 0};
  const Ctds kCtd7 = {CTD_UNUSED, CTD_FCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
      CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_NEXT};
  const std::string kCtdString7 = "[p...].[...]N";

  const Primary kPrimary8 = Primary::FromString("UACGUUGGUGCUUAA");
  const Secondary kSecondary8 = {14, -1, 6, -1, -1, -1, 2, -1, 12, -1, -1, -1, 8, -1, 0};
  const Ctds kCtd8 = {CTD_UNUSED, CTD_NA, CTD_LCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
      CTD_NA, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_LCOAX_WITH_NEXT};
  const std::string kCtdString8 = "[Mp...].[...]mN";

  const Primary kPrimary9 = Primary::FromString("UACGUUGGUGCUUA");
  const Secondary kSecondary9 = {13, -1, 6, -1, -1, -1, 2, -1, 12, -1, -1, -1, 8, 0};
  const Ctds kCtd9 = {CTD_UNUSED, CTD_NA, CTD_RCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
      CTD_NA, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_RCOAX_WITH_NEXT};
  const std::string kCtdString9 = "[mp...]M[...]N";

  const Primary kPrimary10 = Primary::FromString("UACGUUGGUGCUU");
  const Secondary kSecondary10 = {12, 5, -1, -1, -1, 1, -1, 11, -1, -1, -1, 7, 0};
  const Ctds kCtd10 = {CTD_UNUSED, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
      CTD_FCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_PREV};
  const std::string kCtdString10 = "[[...].n...]P";

  const Primary kPrimary11 = Primary::FromString("UACGUUGGUGCUUAA");
  const Secondary kSecondary11 = {14, -1, 6, -1, -1, -1, 2, -1, 12, -1, -1, -1, 8, -1, 0};
  const Ctds kCtd11 = {CTD_UNUSED, CTD_NA, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
      CTD_RCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_RCOAX_WITH_PREV};
  const std::string kCtdString11 = "[M[...].n...]mP";

  const Primary kPrimary12 = Primary::FromString("UACGUUGGUGCUUA");
  const Secondary kSecondary12 = {13, 5, -1, -1, -1, 1, -1, 11, -1, -1, -1, 7, -1, 0};
  const Ctds kCtd12 = {CTD_UNUSED, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
      CTD_LCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_LCOAX_WITH_PREV};
  const std::string kCtdString12 = "[[...]mn...]MP";

  const Primary kPrimary13 = Primary::FromString("UACGUUGGUGCU");
  const Secondary kSecondary13 = {11, 5, -1, -1, -1, 1, 10, -1, -1, -1, 6, 0};
  const Ctds kCtd13 = {CTD_UNUSED, CTD_FCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
      CTD_FCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_UNUSED};
  const std::string kCtdString13 = "[n...]p...]]";

  const Primary kPrimary14 = Primary::FromString("AAAAAAA");
  const Secondary kSecondary14 = {-1, 5, -1, -1, -1, 1, -1};
  const Ctds kCtd14 = {CTD_NA, CTD_MISMATCH, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA};
  const std::string kCtdString14 = "m[...]M";

  const Primary kPrimary15 = Primary::FromString("UACGUUGGUGCUAAAA");
  const Secondary kSecondary15 = {15, -1, 6, -1, -1, -1, 2, -1, -1, 13, -1, -1, -1, 9, -1, 0};
  const Ctds kCtd15 = {CTD_UNUSED, CTD_NA, CTD_MISMATCH, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
      CTD_NA, CTD_MISMATCH, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_UNUSED};
  const std::string kCtdString15 = "[m[...]Mm[...]M]";
};

TEST_F(ParsingTest, SecondaryCtdToCtdString) {
  EXPECT_EQ(kCtdString1, kCtd1.ToString(kSecondary1));
  EXPECT_EQ(kCtdString2, kCtd2.ToString(kSecondary2));
  EXPECT_EQ(kCtdString3, kCtd3.ToString(kSecondary3));
  EXPECT_EQ(kCtdString4, kCtd4.ToString(kSecondary4));
  EXPECT_EQ(kCtdString5, kCtd5.ToString(kSecondary5));
  EXPECT_EQ(kCtdString6, kCtd6.ToString(kSecondary6));
  EXPECT_EQ(kCtdString7, kCtd7.ToString(kSecondary7));
  EXPECT_EQ(kCtdString8, kCtd8.ToString(kSecondary8));
  EXPECT_EQ(kCtdString9, kCtd9.ToString(kSecondary9));
  EXPECT_EQ(kCtdString10, kCtd10.ToString(kSecondary10));
  EXPECT_EQ(kCtdString11, kCtd11.ToString(kSecondary11));
  EXPECT_EQ(kCtdString12, kCtd12.ToString(kSecondary12));
  EXPECT_EQ(kCtdString13, kCtd13.ToString(kSecondary13));
  EXPECT_EQ(kCtdString14, kCtd14.ToString(kSecondary14));
  EXPECT_EQ(kCtdString15, kCtd15.ToString(kSecondary15));
}

TEST_F(ParsingTest, ParseCtd) {
  EXPECT_EQ(std::make_tuple(kPrimary1, kSecondary1, kCtd1),
      ParsePrimaryCtdString(kPrimary1.ToString(), kCtdString1));
  EXPECT_EQ(std::make_tuple(kPrimary2, kSecondary2, kCtd2),
      ParsePrimaryCtdString(kPrimary2.ToString(), kCtdString2));
  EXPECT_EQ(std::make_tuple(kPrimary3, kSecondary3, kCtd3),
      ParsePrimaryCtdString(kPrimary3.ToString(), kCtdString3));
  EXPECT_EQ(std::make_tuple(kPrimary4, kSecondary4, kCtd4),
      ParsePrimaryCtdString(kPrimary4.ToString(), kCtdString4));
  EXPECT_EQ(std::make_tuple(kPrimary5, kSecondary5, kCtd5),
      ParsePrimaryCtdString(kPrimary5.ToString(), kCtdString5));
  EXPECT_EQ(std::make_tuple(kPrimary6, kSecondary6, kCtd6),
      ParsePrimaryCtdString(kPrimary6.ToString(), kCtdString6));
  EXPECT_EQ(std::make_tuple(kPrimary7, kSecondary7, kCtd7),
      ParsePrimaryCtdString(kPrimary7.ToString(), kCtdString7));
  EXPECT_EQ(std::make_tuple(kPrimary8, kSecondary8, kCtd8),
      ParsePrimaryCtdString(kPrimary8.ToString(), kCtdString8));
  EXPECT_EQ(std::make_tuple(kPrimary9, kSecondary9, kCtd9),
      ParsePrimaryCtdString(kPrimary9.ToString(), kCtdString9));
  EXPECT_EQ(std::make_tuple(kPrimary10, kSecondary10, kCtd10),
      ParsePrimaryCtdString(kPrimary10.ToString(), kCtdString10));
  EXPECT_EQ(std::make_tuple(kPrimary11, kSecondary11, kCtd11),
      ParsePrimaryCtdString(kPrimary11.ToString(), kCtdString11));
  EXPECT_EQ(std::make_tuple(kPrimary12, kSecondary12, kCtd12),
      ParsePrimaryCtdString(kPrimary12.ToString(), kCtdString12));
  EXPECT_EQ(std::make_tuple(kPrimary13, kSecondary13, kCtd13),
      ParsePrimaryCtdString(kPrimary13.ToString(), kCtdString13));
  EXPECT_EQ(std::make_tuple(kPrimary14, kSecondary14, kCtd14),
      ParsePrimaryCtdString(kPrimary14.ToString(), kCtdString14));
  EXPECT_EQ(std::make_tuple(kPrimary15, kSecondary15, kCtd15),
      ParsePrimaryCtdString(kPrimary15.ToString(), kCtdString15));
}

TEST_F(ParsingTest, IsCtdString) {
  EXPECT_TRUE(Ctds::IsCtdString(kCtdString1));
  EXPECT_TRUE(Ctds::IsCtdString(kCtdString2));
  EXPECT_TRUE(Ctds::IsCtdString(kCtdString3));
  EXPECT_TRUE(Ctds::IsCtdString(kCtdString4));
  EXPECT_TRUE(Ctds::IsCtdString(kCtdString5));
  EXPECT_TRUE(Ctds::IsCtdString(kCtdString6));
  EXPECT_TRUE(Ctds::IsCtdString(kCtdString7));
  EXPECT_TRUE(Ctds::IsCtdString(kCtdString8));
  EXPECT_TRUE(Ctds::IsCtdString(kCtdString9));
  EXPECT_TRUE(Ctds::IsCtdString(kCtdString10));
  EXPECT_TRUE(Ctds::IsCtdString(kCtdString11));
  EXPECT_TRUE(Ctds::IsCtdString(kCtdString12));
  EXPECT_TRUE(Ctds::IsCtdString(kCtdString13));
  EXPECT_TRUE(Ctds::IsCtdString(kCtdString14));
  EXPECT_TRUE(Ctds::IsCtdString(kCtdString15));
  EXPECT_FALSE(Ctds::IsCtdString("(...)"));
  EXPECT_FALSE(Ctds::IsCtdString("((...))"));
}

}  // namespace mrna
