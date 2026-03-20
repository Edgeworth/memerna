// Copyright 2016 Eliot Courtney.
#ifndef API_MFE_H_
#define API_MFE_H_

#include <variant>

#include "backends/common/base/dp.h"
#include "backends/stack/mfe/dp.h"
#include "model/energy.h"
#include "util/util.h"

namespace mrna::mfe {

// Monostate if there is no DP state supported/used, like with RNAstructure.
using DpState = std::variant<std::monostate, md::base::DpState, md::stack::DpState>;

struct MfeResult {
  DpState dp;
  Energy energy = ZERO_E;
};

[[nodiscard]] inline const md::base::DpState* MaybeGetBaseDpState(const DpState& dp) {
  return std::visit(
      overloaded{
          [](const md::base::DpState& base) -> const md::base::DpState* { return &base; },
          [](const md::stack::DpState& stack) -> const md::base::DpState* { return &stack.base; },
          [](const std::monostate&) -> const md::base::DpState* { return nullptr; },
      },
      dp);
}

}  // namespace mrna::mfe

#endif  // API_MFE_H_
