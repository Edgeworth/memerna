// Copyright 2016 Eliot Courtney.
#ifndef API_MFE_H_
#define API_MFE_H_

#include <variant>

#include "backends/common/base/dp.h"
#include "backends/stack/mfe/mfe.h"
#include "model/energy.h"

namespace mrna::mfe {

// Monostate if there is no DP state supported/used, like with RNAstructure.
using DpState = std::variant<std::monostate, md::base::DpState, md::stack::DpState>;

struct MfeResult {
  DpState dp;
  Energy energy = ZERO_E;
};

}  // namespace mrna::mfe

#endif  // API_MFE_H_
