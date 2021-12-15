// Copyright 2016 E.
#include "model/parsing.h"

#include "gtest/gtest.h"

namespace mrna {

class ParsingTest : public testing::Test {
 public:
  const std::string kCtd1 = "";
  const std::string kCtd2 = ".";
  const std::string kCtd3 = "[[...].[...]3]";
  const std::string kCtd4 = ".n...]mp...]M.";
  const std::string kCtd5 = ".[...]3[...]..";
  const std::string kCtd6 = "n[...]]p...]..";
  const std::string kCtd7 = "[p...].[...]N";
  const std::string kCtd8 = "[Mp...].[...]mN";
  const std::string kCtd9 = "[mp...]M[...]N";
  const std::string kCtd10 = "[[...].n...]P";
  const std::string kCtd11 = "[M[...].n...]mP";
  const std::string kCtd12 = "[[...]mn...]MP";
  const std::string kCtd13 = "[n...]p...]]";
  const std::string kCtd14 = "m[...]M";
  const std::string kCtd15 = "[m[...]Mm[...]M]";

  const computed_t kComputed1 = {{{}, {}}, {}, MAX_E};
  const computed_t kComputed2 = {{{A}, {-1}}, {CTD_NA}, MAX_E};
  const computed_t kComputed3 = {
      {StringToPrimary("UACGUUGGUGCUUA"), {13, 5, -1, -1, -1, 1, -1, 11, -1, -1, -1, 7, -1, 0}},
      {CTD_UNUSED, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_3_DANGLE, CTD_NA, CTD_NA,
          CTD_NA, CTD_NA, CTD_NA, CTD_UNUSED},
      MAX_E};
  const computed_t kComputed4 = {
      {StringToPrimary("UACGUUGGUGCUUA"), {-1, 5, -1, -1, -1, 1, -1, 11, -1, -1, -1, 7, -1, -1}},
      {CTD_NA, CTD_RCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_RCOAX_WITH_PREV,
          CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA},
      MAX_E};
  const computed_t kComputed5 = {
      {StringToPrimary("UACGUUGGUGCUUA"), {-1, 5, -1, -1, -1, 1, -1, 11, -1, -1, -1, 7, -1, -1}},
      {CTD_NA, CTD_3_DANGLE, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_UNUSED, CTD_NA, CTD_NA,
          CTD_NA, CTD_NA, CTD_NA, CTD_NA},
      MAX_E};
  const computed_t kComputed6 = {
      {StringToPrimary("UACGUUGGUGCUUA"), {6, 5, -1, -1, -1, 1, 0, 11, -1, -1, -1, 7, -1, -1}},
      {CTD_FCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_PREV,
          CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA},
      MAX_E};
  const computed_t kComputed7 = {
      {StringToPrimary("UACGUUGGUGCUU"), {12, 5, -1, -1, -1, 1, -1, 11, -1, -1, -1, 7, 0}},
      {CTD_UNUSED, CTD_FCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_UNUSED, CTD_NA,
          CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_NEXT},
      MAX_E};
  const computed_t kComputed8 = {{StringToPrimary("UACGUUGGUGCUUAA"),
                                     {14, -1, 6, -1, -1, -1, 2, -1, 12, -1, -1, -1, 8, -1, 0}},
      {CTD_UNUSED, CTD_NA, CTD_LCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_UNUSED,
          CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_LCOAX_WITH_NEXT},
      MAX_E};
  const computed_t kComputed9 = {
      {StringToPrimary("UACGUUGGUGCUUA"), {13, -1, 6, -1, -1, -1, 2, -1, 12, -1, -1, -1, 8, 0}},
      {CTD_UNUSED, CTD_NA, CTD_RCOAX_WITH_PREV, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_UNUSED,
          CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_RCOAX_WITH_NEXT},
      MAX_E};
  const computed_t kComputed10 = {
      {StringToPrimary("UACGUUGGUGCUU"), {12, 5, -1, -1, -1, 1, -1, 11, -1, -1, -1, 7, 0}},
      {CTD_UNUSED, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_NEXT, CTD_NA,
          CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_PREV},
      MAX_E};
  const computed_t kComputed11 = {{StringToPrimary("UACGUUGGUGCUUAA"),
                                      {14, -1, 6, -1, -1, -1, 2, -1, 12, -1, -1, -1, 8, -1, 0}},
      {CTD_UNUSED, CTD_NA, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_RCOAX_WITH_NEXT,
          CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_RCOAX_WITH_PREV},
      MAX_E};
  const computed_t kComputed12 = {
      {StringToPrimary("UACGUUGGUGCUUA"), {13, 5, -1, -1, -1, 1, -1, 11, -1, -1, -1, 7, -1, 0}},
      {CTD_UNUSED, CTD_UNUSED, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_LCOAX_WITH_NEXT, CTD_NA,
          CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_LCOAX_WITH_PREV},
      MAX_E};
  const computed_t kComputed13 = {
      {StringToPrimary("UACGUUGGUGCU"), {11, 5, -1, -1, -1, 1, 10, -1, -1, -1, 6, 0}},
      {CTD_UNUSED, CTD_FCOAX_WITH_NEXT, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_FCOAX_WITH_PREV, CTD_NA,
          CTD_NA, CTD_NA, CTD_NA, CTD_UNUSED},
      MAX_E};
  const computed_t kComputed14 = {{StringToPrimary("AAAAAAA"), {-1, 5, -1, -1, -1, 1, -1}},
      {CTD_NA, CTD_MISMATCH, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA}, MAX_E};
  const computed_t kComputed15 = {{StringToPrimary("UACGUUGGUGCUAAAA"),
                                      {15, -1, 6, -1, -1, -1, 2, -1, -1, 13, -1, -1, -1, 9, -1, 0}},
      {CTD_UNUSED, CTD_NA, CTD_MISMATCH, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA,
          CTD_MISMATCH, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_NA, CTD_UNUSED},
      MAX_E};
};

TEST_F(ParsingTest, ComputedToCtdString) {
  EXPECT_EQ(kCtd1, ComputedToCtdString(kComputed1));
  EXPECT_EQ(kCtd2, ComputedToCtdString(kComputed2));
  EXPECT_EQ(kCtd3, ComputedToCtdString(kComputed3));
  EXPECT_EQ(kCtd4, ComputedToCtdString(kComputed4));
  EXPECT_EQ(kCtd5, ComputedToCtdString(kComputed5));
  EXPECT_EQ(kCtd6, ComputedToCtdString(kComputed6));
  EXPECT_EQ(kCtd7, ComputedToCtdString(kComputed7));
  EXPECT_EQ(kCtd8, ComputedToCtdString(kComputed8));
  EXPECT_EQ(kCtd9, ComputedToCtdString(kComputed9));
  EXPECT_EQ(kCtd10, ComputedToCtdString(kComputed10));
  EXPECT_EQ(kCtd11, ComputedToCtdString(kComputed11));
  EXPECT_EQ(kCtd12, ComputedToCtdString(kComputed12));
  EXPECT_EQ(kCtd13, ComputedToCtdString(kComputed13));
  EXPECT_EQ(kCtd14, ComputedToCtdString(kComputed14));
  EXPECT_EQ(kCtd15, ComputedToCtdString(kComputed15));
}

TEST_F(ParsingTest, ParseCtdComputed) {
  EXPECT_EQ(kComputed1, ParseCtdComputed(PrimaryToString(kComputed1.s.r), kCtd1));
  EXPECT_EQ(kComputed2, ParseCtdComputed(PrimaryToString(kComputed2.s.r), kCtd2));
  EXPECT_EQ(kComputed3, ParseCtdComputed(PrimaryToString(kComputed3.s.r), kCtd3));
  EXPECT_EQ(kComputed4, ParseCtdComputed(PrimaryToString(kComputed4.s.r), kCtd4));
  EXPECT_EQ(kComputed5, ParseCtdComputed(PrimaryToString(kComputed5.s.r), kCtd5));
  EXPECT_EQ(kComputed6, ParseCtdComputed(PrimaryToString(kComputed6.s.r), kCtd6));
  EXPECT_EQ(kComputed7, ParseCtdComputed(PrimaryToString(kComputed7.s.r), kCtd7));
  EXPECT_EQ(kComputed8, ParseCtdComputed(PrimaryToString(kComputed8.s.r), kCtd8));
  EXPECT_EQ(kComputed9, ParseCtdComputed(PrimaryToString(kComputed9.s.r), kCtd9));
  EXPECT_EQ(kComputed10, ParseCtdComputed(PrimaryToString(kComputed10.s.r), kCtd10));
  EXPECT_EQ(kComputed11, ParseCtdComputed(PrimaryToString(kComputed11.s.r), kCtd11));
  EXPECT_EQ(kComputed12, ParseCtdComputed(PrimaryToString(kComputed12.s.r), kCtd12));
  EXPECT_EQ(kComputed13, ParseCtdComputed(PrimaryToString(kComputed13.s.r), kCtd13));
  EXPECT_EQ(kComputed14, ParseCtdComputed(PrimaryToString(kComputed14.s.r), kCtd14));
  EXPECT_EQ(kComputed15, ParseCtdComputed(PrimaryToString(kComputed15.s.r), kCtd15));
}

TEST_F(ParsingTest, IsCtdString) {
  EXPECT_TRUE(IsCtdString(kCtd1));
  EXPECT_TRUE(IsCtdString(kCtd2));
  EXPECT_TRUE(IsCtdString(kCtd3));
  EXPECT_TRUE(IsCtdString(kCtd4));
  EXPECT_TRUE(IsCtdString(kCtd5));
  EXPECT_TRUE(IsCtdString(kCtd6));
  EXPECT_TRUE(IsCtdString(kCtd7));
  EXPECT_TRUE(IsCtdString(kCtd8));
  EXPECT_TRUE(IsCtdString(kCtd9));
  EXPECT_TRUE(IsCtdString(kCtd10));
  EXPECT_TRUE(IsCtdString(kCtd11));
  EXPECT_TRUE(IsCtdString(kCtd12));
  EXPECT_TRUE(IsCtdString(kCtd13));
  EXPECT_TRUE(IsCtdString(kCtd14));
  EXPECT_TRUE(IsCtdString(kCtd15));
  EXPECT_FALSE(IsCtdString("(...)"));
  EXPECT_FALSE(IsCtdString("((...))"));
}

}  // namespace mrna
