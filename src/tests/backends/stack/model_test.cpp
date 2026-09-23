// Copyright 2016 Eliot Courtney.
#include "backends/stack/energy/model.h"

#include "api/ctx/backend_cfg.h"
#include "backends/stack/energy/boltz_model.h"
#include "gtest/gtest.h"
#include "model/energy.h"
#include "tests/init.h"

namespace mrna::md::stack {

class ModelTestStack : public testing::Test {};

TEST_F(ModelTestStack, IsValid) {
  for (const auto& m : stack_ms) EXPECT_TRUE(m->IsValid());

  auto random_model = Model::Random(kT04Cfg, RandomModelCfg(0, E(0.1), E(10.0)));
  EXPECT_TRUE(random_model->IsValid());
  random_model->multiloop_c = E(1.0);
  EXPECT_TRUE(random_model->IsValid());
}

TEST_F(ModelTestStack, RandomAndBoltzModelsFillPenultimateStack) {
  const auto model = Model::Random(kT22Cfg, RandomModelCfg(0, E(1.0), E(1.0)));
  const auto boltz_model = BoltzModel::Create(model);

  for (int a = 0; a < 4; ++a) {
    for (int b = 0; b < 4; ++b) {
      for (int c = 0; c < 4; ++c) {
        for (int d = 0; d < 4; ++d) {
          EXPECT_EQ(E(1.0), model->penultimate_stack[a][b][c][d]);
          EXPECT_EQ(Boltz(E(1.0)), boltz_model->penultimate_stack[a][b][c][d]);
        }
      }
    }
  }
}

}  // namespace mrna::md::stack
