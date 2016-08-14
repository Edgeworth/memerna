#include "fold.h"
#include "parsing.h"
#include "gtest/gtest.h"

namespace memerna {
namespace fold {

class FoldTest : public testing::Test {
public:
  rna_t kStacking = parsing::StringToRna("GGGGAAACCCC");
  rna_t kTest = parsing::StringToRna("UUGAAAAGCGGUUCCGUUCAGUCCUACUCACACGUCCGUCACACAUUAUGCCGGUAGAUA");
};


TEST_F(FoldTest, T04) {
  EXPECT_EQ(-45, Fold(kStacking));
}

}
}
