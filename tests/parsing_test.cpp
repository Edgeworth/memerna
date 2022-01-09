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

  const Primary kPrimary2 = StringToPrimary("A");
  const Secondary kSecondary2 = {-1};
  const Ctds kCtd2 = {CTD_NA};
  const std::string kCtdString2 = ".";

  const Primary kPrimary3 = StringToPrimary("UACGUUGGUGCUUA");
  const Secondary kSecondary3 = {13, 5, -1, -1, -1, 1, -1, 11, -1, -1, -1, 7, -1, 0};
  const Ctds kCtd3 = {CTD_UNUSED, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_3_DANGLE,
      CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_UNUSED};
  const std::string kCtdString3 = "[[...].[...]3]";

  const Primary kPrimary4 = StringToPrimary("UACGUUGGUGCUUA");
  const Secondary kSecondary4 = {-1, 5, -1, -1, -1, 1, -1, 11, -1, -1, -1, 7, -1, -1};
  const Ctds kCtd4 = {CTD_NA, CTD_RCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
      CTD_RCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA};
  const std::string kCtdString4 = ".n...]mp...]M.";

  const Primary kPrimary5 = StringToPrimary("UACGUUGGUGCUUA");
  const Secondary kSecondary5 = {-1, 5, -1, -1, -1, 1, -1, 11, -1, -1, -1, 7, -1, -1};
  const Ctds kCtd5 = {CTD_NA, CTD_3_DANGLE, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_UNUSED,
      CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA};
  const std::string kCtdString5 = ".[...]3[...]..";

  const Primary kPrimary6 = StringToPrimary("UACGUUGGUGCUUA");
  const Secondary kSecondary6 = {6, 5, -1, -1, -1, 1, 0, 11, -1, -1, -1, 7, -1, -1};
  const Ctds kCtd6 = {CTD_FCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
      CTD_FCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA};
  const std::string kCtdString6 = "n[...]]p...]..";

  const Primary kPrimary7 = StringToPrimary("UACGUUGGUGCUU");
  const Secondary kSecondary7 = {12, 5, -1, -1, -1, 1, -1, 11, -1, -1, -1, 7, 0};
  const Ctds kCtd7 = {CTD_UNUSED, CTD_FCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
      CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_NEXT};
  const std::string kCtdString7 = "[p...].[...]N";

  const Primary kPrimary8 = StringToPrimary("UACGUUGGUGCUUAA");
  const Secondary kSecondary8 = {14, -1, 6, -1, -1, -1, 2, -1, 12, -1, -1, -1, 8, -1, 0};
  const Ctds kCtd8 = {CTD_UNUSED, CTD_NA, CTD_LCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
      CTD_NA, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_LCOAX_WITH_NEXT};
  const std::string kCtdString8 = "[Mp...].[...]mN";

  const Primary kPrimary9 = StringToPrimary("UACGUUGGUGCUUA");
  const Secondary kSecondary9 = {13, -1, 6, -1, -1, -1, 2, -1, 12, -1, -1, -1, 8, 0};
  const Ctds kCtd9 = {CTD_UNUSED, CTD_NA, CTD_RCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
      CTD_NA, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_RCOAX_WITH_NEXT};
  const std::string kCtdString9 = "[mp...]M[...]N";

  const Primary kPrimary10 = StringToPrimary("UACGUUGGUGCUU");
  const Secondary kSecondary10 = {12, 5, -1, -1, -1, 1, -1, 11, -1, -1, -1, 7, 0};
  const Ctds kCtd10 = {CTD_UNUSED, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
      CTD_FCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_PREV};
  const std::string kCtdString10 = "[[...].n...]P";

  const Primary kPrimary11 = StringToPrimary("UACGUUGGUGCUUAA");
  const Secondary kSecondary11 = {14, -1, 6, -1, -1, -1, 2, -1, 12, -1, -1, -1, 8, -1, 0};
  const Ctds kCtd11 = {CTD_UNUSED, CTD_NA, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
      CTD_RCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_RCOAX_WITH_PREV};
  const std::string kCtdString11 = "[M[...].n...]mP";

  const Primary kPrimary12 = StringToPrimary("UACGUUGGUGCUUA");
  const Secondary kSecondary12 = {13, 5, -1, -1, -1, 1, -1, 11, -1, -1, -1, 7, -1, 0};
  const Ctds kCtd12 = {CTD_UNUSED, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
      CTD_LCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_LCOAX_WITH_PREV};
  const std::string kCtdString12 = "[[...]mn...]MP";

  const Primary kPrimary13 = StringToPrimary("UACGUUGGUGCU");
  const Secondary kSecondary13 = {11, 5, -1, -1, -1, 1, 10, -1, -1, -1, 6, 0};
  const Ctds kCtd13 = {CTD_UNUSED, CTD_FCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
      CTD_FCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_UNUSED};
  const std::string kCtdString13 = "[n...]p...]]";

  const Primary kPrimary14 = StringToPrimary("AAAAAAA");
  const Secondary kSecondary14 = {-1, 5, -1, -1, -1, 1, -1};
  const Ctds kCtd14 = {CTD_NA, CTD_MISMATCH, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA};
  const std::string kCtdString14 = "m[...]M";

  const Primary kPrimary15 = StringToPrimary("UACGUUGGUGCUAAAA");
  const Secondary kSecondary15 = {15, -1, 6, -1, -1, -1, 2, -1, -1, 13, -1, -1, -1, 9, -1, 0};
  const Ctds kCtd15 = {CTD_UNUSED, CTD_NA, CTD_MISMATCH, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
      CTD_NA, CTD_MISMATCH, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_UNUSED};
  const std::string kCtdString15 = "[m[...]Mm[...]M]";
};

TEST_F(ParsingTest, ComputedToCtdString) {
  EXPECT_EQ(kCtdString1, CtdString(kSecondary1, kCtd1));
  EXPECT_EQ(kCtdString2, CtdString(kSecondary2, kCtd2));
  EXPECT_EQ(kCtdString3, CtdString(kSecondary3, kCtd3));
  EXPECT_EQ(kCtdString4, CtdString(kSecondary4, kCtd4));
  EXPECT_EQ(kCtdString5, CtdString(kSecondary5, kCtd5));
  EXPECT_EQ(kCtdString6, CtdString(kSecondary6, kCtd6));
  EXPECT_EQ(kCtdString7, CtdString(kSecondary7, kCtd7));
  EXPECT_EQ(kCtdString8, CtdString(kSecondary8, kCtd8));
  EXPECT_EQ(kCtdString9, CtdString(kSecondary9, kCtd9));
  EXPECT_EQ(kCtdString10, CtdString(kSecondary10, kCtd10));
  EXPECT_EQ(kCtdString11, CtdString(kSecondary11, kCtd11));
  EXPECT_EQ(kCtdString12, CtdString(kSecondary12, kCtd12));
  EXPECT_EQ(kCtdString13, CtdString(kSecondary13, kCtd13));
  EXPECT_EQ(kCtdString14, CtdString(kSecondary14, kCtd14));
  EXPECT_EQ(kCtdString15, CtdString(kSecondary15, kCtd15));
}

TEST_F(ParsingTest, ParseCtdComputed) {
  EXPECT_EQ(std::make_tuple(kPrimary1, kSecondary1, kCtd1),
      ParsePrimaryCtdString(PrimaryToString(kPrimary1), kCtdString1));
  EXPECT_EQ(std::make_tuple(kPrimary2, kSecondary2, kCtd2),
      ParsePrimaryCtdString(PrimaryToString(kPrimary2), kCtdString2));
  EXPECT_EQ(std::make_tuple(kPrimary3, kSecondary3, kCtd3),
      ParsePrimaryCtdString(PrimaryToString(kPrimary3), kCtdString3));
  EXPECT_EQ(std::make_tuple(kPrimary4, kSecondary4, kCtd4),
      ParsePrimaryCtdString(PrimaryToString(kPrimary4), kCtdString4));
  EXPECT_EQ(std::make_tuple(kPrimary5, kSecondary5, kCtd5),
      ParsePrimaryCtdString(PrimaryToString(kPrimary5), kCtdString5));
  EXPECT_EQ(std::make_tuple(kPrimary6, kSecondary6, kCtd6),
      ParsePrimaryCtdString(PrimaryToString(kPrimary6), kCtdString6));
  EXPECT_EQ(std::make_tuple(kPrimary7, kSecondary7, kCtd7),
      ParsePrimaryCtdString(PrimaryToString(kPrimary7), kCtdString7));
  EXPECT_EQ(std::make_tuple(kPrimary8, kSecondary8, kCtd8),
      ParsePrimaryCtdString(PrimaryToString(kPrimary8), kCtdString8));
  EXPECT_EQ(std::make_tuple(kPrimary9, kSecondary9, kCtd9),
      ParsePrimaryCtdString(PrimaryToString(kPrimary9), kCtdString9));
  EXPECT_EQ(std::make_tuple(kPrimary10, kSecondary10, kCtd10),
      ParsePrimaryCtdString(PrimaryToString(kPrimary10), kCtdString10));
  EXPECT_EQ(std::make_tuple(kPrimary11, kSecondary11, kCtd11),
      ParsePrimaryCtdString(PrimaryToString(kPrimary11), kCtdString11));
  EXPECT_EQ(std::make_tuple(kPrimary12, kSecondary12, kCtd12),
      ParsePrimaryCtdString(PrimaryToString(kPrimary12), kCtdString12));
  EXPECT_EQ(std::make_tuple(kPrimary13, kSecondary13, kCtd13),
      ParsePrimaryCtdString(PrimaryToString(kPrimary13), kCtdString13));
  EXPECT_EQ(std::make_tuple(kPrimary14, kSecondary14, kCtd14),
      ParsePrimaryCtdString(PrimaryToString(kPrimary14), kCtdString14));
  EXPECT_EQ(std::make_tuple(kPrimary15, kSecondary15, kCtd15),
      ParsePrimaryCtdString(PrimaryToString(kPrimary15), kCtdString15));
}

TEST_F(ParsingTest, IsCtdString) {
  EXPECT_TRUE(IsCtdString(kCtdString1));
  EXPECT_TRUE(IsCtdString(kCtdString2));
  EXPECT_TRUE(IsCtdString(kCtdString3));
  EXPECT_TRUE(IsCtdString(kCtdString4));
  EXPECT_TRUE(IsCtdString(kCtdString5));
  EXPECT_TRUE(IsCtdString(kCtdString6));
  EXPECT_TRUE(IsCtdString(kCtdString7));
  EXPECT_TRUE(IsCtdString(kCtdString8));
  EXPECT_TRUE(IsCtdString(kCtdString9));
  EXPECT_TRUE(IsCtdString(kCtdString10));
  EXPECT_TRUE(IsCtdString(kCtdString11));
  EXPECT_TRUE(IsCtdString(kCtdString12));
  EXPECT_TRUE(IsCtdString(kCtdString13));
  EXPECT_TRUE(IsCtdString(kCtdString14));
  EXPECT_TRUE(IsCtdString(kCtdString15));
  EXPECT_FALSE(IsCtdString("(...)"));
  EXPECT_FALSE(IsCtdString("((...))"));
}

}  // namespace mrna
