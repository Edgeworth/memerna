#ifndef MODEL_SECONDARY_H_
#define MODEL_SECONDARY_H_

#include "model/primary.h"

namespace mrna {

struct Secondary {
  Secondary() = default;
  Secondary(const Secondary&) = default;
  Secondary(Secondary&&) = default;
  Secondary& operator=(Secondary&&) = default;

  explicit Secondary(const Primary& r) : r(r), p(r.size(), -1) {}
  Secondary(const Primary& r, const std::vector<int>& p) : r(r), p(p) {}

  bool operator==(const Secondary& o) const { return r == o.r && p == o.p; }
  bool operator!=(const Secondary& o) const { return !(*this == o); }
  bool operator<(const Secondary& o) const { return std::tie(r, p) < std::tie(o.r, o.p); }

  Primary r;
  std::vector<int> p;
};

}  // namespace mrna

#endif  // MODEL_SECONDARY_H_
