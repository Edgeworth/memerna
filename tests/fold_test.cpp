#include <cstdlib>
#include "fold/fold.h"
#include "fold/fold_globals.h"
#include "parsing.h"
#include "gtest/gtest.h"
#include "common_test.h"

namespace memerna {
namespace fold {

class FoldTest : public testing::Test {
public:
  rna_t kStacking = parsing::StringToRna("GGGGAAACCCC");
  rna_t kTest = parsing::StringToRna("UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA");
};


TEST_F(FoldTest, T04) {
  if (EnergyModelChecksum() != T04_MODEL_HASH) {
    printf("Skipping energy model specific energy tests.");
    return;
  }

  EXPECT_EQ(-45, Fold(kStacking));
  EXPECT_EQ(-51, Fold(kTest));
}

TEST_F(FoldTest, Constants) {
  if (EnergyModelChecksum() != T04_MODEL_HASH) {
    printf("Skipping energy model specific energy tests.");
    return;
  }

  r = kStacking;
  fold::InitFold();
  EXPECT_EQ(-21 - 4 - 16, g_min_mismatch_coax);
  EXPECT_EQ(-34, g_min_flush_coax);
  EXPECT_EQ(-26, g_min_twoloop_not_stack);

  energy_t augubranch[4][4] = {
      {-6, -6, -6, 5 - 6},
      {-6, -6, -6, -6},
      {-6, -6, -6, 5 - 6},
      {5 - 6, -6, 5 - 6, -6}
  };
  EXPECT_EQ(sizeof(augubranch), sizeof(g_augubranch));
  EXPECT_TRUE(std::memcmp(augubranch, g_augubranch, sizeof(augubranch)) == 0);
}

TEST_F(FoldTest, Helpers) {
  EXPECT_EQ(0, fold::MaxNumContiguous(parsing::StringToRna("")));
  EXPECT_EQ(1, fold::MaxNumContiguous(parsing::StringToRna("A")));
  EXPECT_EQ(2, fold::MaxNumContiguous(parsing::StringToRna("AA")));
  EXPECT_EQ(2, fold::MaxNumContiguous(parsing::StringToRna("GUAAC")));
  EXPECT_EQ(1, fold::MaxNumContiguous(parsing::StringToRna("GUACA")));
  EXPECT_EQ(3, fold::MaxNumContiguous(parsing::StringToRna("GAUCCC")));
  EXPECT_EQ(3, fold::MaxNumContiguous(parsing::StringToRna("GGGAUC")));
  EXPECT_EQ(4, fold::MaxNumContiguous(parsing::StringToRna("GGGAUCAAAA")));
  EXPECT_EQ(5, fold::MaxNumContiguous(parsing::StringToRna("GGGAUUUUUCAAAA")));
}

}
}
