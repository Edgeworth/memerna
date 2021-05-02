// Copyright 2016, E.
//
// This file is part of memerna.
//
// memerna is free software: you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// memerna is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
// the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with memerna.
// If not, see <http://www.gnu.org/licenses/>.
#include <cstdlib>

#include "common_test.h"
#include "context.h"
#include "parsing.h"

namespace memerna {
namespace fold {

class FoldAlgTest : public testing::TestWithParam<context_opt_t::TableAlg> {
public:
  energy_t FoldEnergy(const std::string& s) {
    return Context(parsing::StringToPrimary(s), g_em, context_opt_t(GetParam())).Fold().energy;
  }
};

TEST_P(FoldAlgTest, T04) {
  ONLY_FOR_THIS_MODEL(g_em, T04_MODEL_HASH);

  EXPECT_EQ(-45, FoldEnergy("GGGGAAACCCC"));
  EXPECT_EQ(-51, FoldEnergy("UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA"));
  EXPECT_EQ(
      -133, FoldEnergy("AAAAACUAGCAUCUAGUGGGCUCCCGAUCGCCUCCUUCUCGUAUUACGUUAAUGCAACUCAAGUGAGCCCGU"));
  EXPECT_EQ(-57, FoldEnergy("AAAUUCCGCUUGACAGCUCGCCACAACGGCAGGAC"));
  EXPECT_EQ(-121, FoldEnergy("ACCGUCAGCUACCGCCGACUAUACUCUUUAGUCAGACGGGG"));
  EXPECT_EQ(-74, FoldEnergy("CACACUCCCGCAAAUGCCGAGUAUCAGAUUACUCCCCGGG"));
  EXPECT_EQ(-32, FoldEnergy("CCCAACGGAGUAACUUAGCGAAUAGCAGGGG"));
  EXPECT_EQ(-40, FoldEnergy("CCGAAGGGGCUGCGGCG"));
  EXPECT_EQ(-120, FoldEnergy("CCGGGCCAGCCCGCUCCUACGGGGGGUC"));
  EXPECT_EQ(-74, FoldEnergy("CCGUAAAGUCGAACCAGACGUGCAUGAGCAAGCGG"));
  EXPECT_EQ(-6, FoldEnergy("CCUCCGGG"));
  EXPECT_EQ(-30, FoldEnergy("CCUGGAUAUUCCGAUGAGCACGUGCGAGGGC"));
  EXPECT_EQ(-65, FoldEnergy("CGCAGGGUCGGACCCGGGAGAACCGCGA"));
  EXPECT_EQ(-60, FoldEnergy("CGCUUAAGGCUAUUUGGCCGGAUCUCCAAGGCA"));
  EXPECT_EQ(-6, FoldEnergy("CGGAAACGG"));
  EXPECT_EQ(-22, FoldEnergy("CUGAAACUGGAAACAGAAAUG"));
  EXPECT_EQ(-12, FoldEnergy("CUUAUAGUUAAGG"));
  EXPECT_EQ(-122, FoldEnergy("GAUGAGGGGAAACGGUGACUGGGACUCAGACAACGAUAGCAGCCAAAUAGGGAAGCUUCCUUC"));
  EXPECT_EQ(-29, FoldEnergy("GCCAAGGCCCCACCCGGA"));
  EXPECT_EQ(-39, FoldEnergy("GCCAGACAAACACGAUUCUUUGAUAGUACUGACUAUUCUACAAUUAGGCC"));
  EXPECT_EQ(-67, FoldEnergy("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCC"));
  EXPECT_EQ(-276,
      FoldEnergy("GCGACCGGGGCUGGCUUGGUAAUGGUACUCCCCUGUCACGGGAGAGAAUGUGGGUUCAAAUCCCAUCGGUCGCGCCA"));
  EXPECT_EQ(-53, FoldEnergy("GCGCCCCAGUCGACGCUGAGCUCCUCUGCU"));
  EXPECT_EQ(-157, FoldEnergy("GGCACAUACUGGGACAACAUUCGUUGGGUUCCCGGGUCGAACGGCAGCCG"));
  EXPECT_EQ(-49, FoldEnergy("GGCCGAUGGCAGCGAUAGC"));
  EXPECT_EQ(-44, FoldEnergy("GGCGCACGCGUUAGCCGGGGAUCCACAGUGC"));
  EXPECT_EQ(-45, FoldEnergy("GGGGAAACCCC"));
  EXPECT_EQ(-29, FoldEnergy("GUCAUGCACUACUGCGAUUCAUACGGAAACAGACG"));
  EXPECT_EQ(-23, FoldEnergy("UACCCUGUUCAGCAUUGGAAAUUUCCUGGG"));
  EXPECT_EQ(-80, FoldEnergy("UCCACGGCUCGACGGCGCACUUAGUGCGUGGG"));
  EXPECT_EQ(-4, FoldEnergy("UGCAAAGCAA"));
  EXPECT_EQ(-208, FoldEnergy("UGGGGAAGUGCCGAUGCGGUACUAUUAUCCACUGUCUAUGGAUAAGUCCCCCGACCU"));
}

INSTANTIATE_TEST_CASE_P(FoldAlgTest, FoldAlgTest, testing::ValuesIn(context_opt_t::TABLE_ALGS));
}  // namespace fold
}  // namespace memerna
