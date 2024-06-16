// Copyright 2016 Eliot Courtney.
#ifndef API_PFN_H_
#define API_PFN_H_

#include <variant>

#include "backends/common/base/dp.h"
#include "model/pfn.h"

namespace mrna::pfn {

using PfnState = std::variant<std::monostate, md::base::PfnState>;

struct PfnResult {
  PfnState state;
  PfnTables pfn;
};

}  // namespace mrna::pfn

#endif  // API_PFN_H_
