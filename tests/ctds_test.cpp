#include "constants.h"
#include "parsing.h"
#include "gtest/gtest.h"
#include "energy/energy_internal.h"

namespace memerna {
namespace energy {

struct ctd_test_t {
  secondary_t secondary;
  internal::branch_ctd_t branch_ctds;
  std::vector<Ctd> base_ctds;
  std::deque<int> branches;
} CTD_TESTS[] = {
    {{{}, {}}, {}, {}, {}}
};

class CtdsTest : public testing::TestWithParam<ctd_test_t> {
};

TEST_P(CtdsTest, ReadWriteRead) {

}

TEST_P(CtdsTest, WriteReadWrite) {
}

INSTANTIATE_TEST_CASE_P(CtdsTest, CtdsTest, testing::ValuesIn(CTD_TESTS));

}
}
