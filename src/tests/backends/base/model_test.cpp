// Copyright 2016 Eliot Courtney.
#include "gtest/gtest.h"
#include "tests/init.h"

namespace mrna::md::base {

class ModelTestBase : public testing::Test {};

TEST_F(ModelTestBase, IsValid) {
  for (const auto& m : base_ms) EXPECT_TRUE(m->IsValid());

  auto random_model = Model::Random(0);
  EXPECT_TRUE(random_model->IsValid());
  random_model->multiloop_c = E(1.0);
  EXPECT_TRUE(random_model->IsValid());
}

}  // namespace mrna::md::base
