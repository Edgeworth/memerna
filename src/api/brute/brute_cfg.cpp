// Copyright 2022 Eliot Courtney.
#include "api/brute/brute_cfg.h"

#include <ostream>

namespace mrna::brute {

std::ostream& operator<<(std::ostream& str, const BruteCfg& o) {
  return str << "BruteCfg{mfe=" << o.mfe << ", subopt=" << o.subopt << ", pfn=" << o.pfn
             << ", subopt_cfg=" << o.subopt_cfg << "}";
}

}  // namespace mrna::brute
