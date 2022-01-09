#ifndef MODEL_PRIMARY_H_
#define MODEL_PRIMARY_H_

#include "model/base.h"

namespace mrna {

using Primary = std::vector<Base>;

Primary GenerateRandomPrimary(int length);
Primary StringToPrimary(const std::string& s);
std::string PrimaryToString(const Primary& r);

}  // namespace mrna

#endif  // MODEL_PRIMARY_H_
