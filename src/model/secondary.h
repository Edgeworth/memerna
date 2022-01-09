#ifndef MODEL_SECONDARY_H_
#define MODEL_SECONDARY_H_

#include "model/primary.h"

namespace mrna {

struct Secondary {
  Secondary() = default;
  Secondary(const Secondary&) = default;
  Secondary(Secondary&&) = default;
  Secondary& operator=(Secondary&&) = default;

  explicit Secondary(Primary r_) : r(std::move(r_)), p(r.size(), -1) {}
  Secondary(Primary r_, std::vector<int> p_) : r(std::move(r_)), p(std::move(p_)) {}

  bool operator==(const Secondary& o) const { return r == o.r && p == o.p; }
  bool operator!=(const Secondary& o) const { return !(*this == o); }
  bool operator<(const Secondary& o) const { return std::tie(r, p) < std::tie(o.r, o.p); }

  Primary r;
  std::vector<int> p;
};

Secondary ParseDotBracketSecondary(const std::string& prim_str, const std::string& pairs_str);
std::vector<int> DotBracketToPairs(const std::string& pairs_str);
std::string PairsToDotBracket(const std::vector<int>& pairs);

}  // namespace mrna

#endif  // MODEL_SECONDARY_H_
