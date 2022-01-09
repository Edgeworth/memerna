#ifndef MODEL_SECONDARY_H_
#define MODEL_SECONDARY_H_

#include "model/primary.h"

namespace mrna {

using Secondary = std::vector<int>;

std::tuple<Primary, Secondary> ParsePrimaryDotBracket(
    const std::string& prim_str, const std::string& pairs_str);
Secondary DotBracketToSecondary(const std::string& pairs_str);
std::string SecondaryToDotBracket(const Secondary& s);

}  // namespace mrna

#endif  // MODEL_SECONDARY_H_
