// Copyright 2016 E.
#ifndef API_MFE_H_
#define API_MFE_H_

#include <variant>

#include "model/energy.h"
#include "models/t04/mfe/dp.h"
#include "models/t22/mfe/mfe.h"

namespace mrna::mfe {

// Monostate if there is no DP state supported/used, like with RNAstructure.
using DpState = std::variant<std::monostate, md::t04::DpState, md::t22::DpState>;

struct MfeResult {
  DpState dp;
  Energy energy = ZERO_E;
};

}  // namespace mrna::mfe

#endif  // API_MFE_H_
