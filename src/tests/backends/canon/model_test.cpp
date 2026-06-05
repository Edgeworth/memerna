// Copyright 2026 Eliot Courtney.
#include "gtest/gtest.h"
#include "tests/init.h"

namespace mrna::md::canon {

class ModelTestCanon : public testing::Test {};

TEST_F(ModelTestCanon, IsValid) {
  for (const auto& m : canon_ms) EXPECT_TRUE(m->IsValid());

  auto random_model = Model::Random(kT04Cfg, RandomModelCfg(0, E(-10.0), E(-0.1)));
  EXPECT_TRUE(random_model->IsValid());
  random_model->multiloop_c = E(1.0);
  EXPECT_TRUE(random_model->IsValid());
}

}  // namespace mrna::md::canon
