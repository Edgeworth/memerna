// Copyright 2026 Eliot Courtney.
#include "api/ctx/backend.h"
#include "api/energy/energy_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "gtest/gtest.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "tests/init.h"

namespace mrna::md::canon {

class EnergyTestCanon : public testing::Test {};

TEST_F(EnergyTestCanon, EquivalentToBase) {
  EXPECT_TRUE(IsEquivalent(*base_t04, *canon_t04));

  const auto cases = {
      kNNDBHairpin1,
      kNNDBBulge1,
      kNNDBInternal2x3,
      kInternal1,
  };

  for (const auto& c : cases) {
    const auto& r = std::get<Primary>(c);
    const auto& s = std::get<Secondary>(c);
    EXPECT_EQ(base_t04->TotalEnergy<erg::EnergyCfg{}>(r, s, nullptr, erg::PseudofreeCfg{}).energy,
        canon_t04->TotalEnergy<erg::EnergyCfg{}>(r, s, nullptr, erg::PseudofreeCfg{}).energy);
    EXPECT_EQ(base_t04
                  ->SubEnergy<erg::EnergyCfg{}>(
                      r, s, nullptr, erg::PseudofreeCfg{}, 0, static_cast<int>(r.size()) - 1)
                  .energy,
        canon_t04
            ->SubEnergy<erg::EnergyCfg{}>(
                r, s, nullptr, erg::PseudofreeCfg{}, 0, static_cast<int>(r.size()) - 1)
            .energy);
  }
}

}  // namespace mrna::md::canon
