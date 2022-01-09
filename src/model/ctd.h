#ifndef MODEL_CTD_H_
#define MODEL_CTD_H_

#include "model/model.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/macros.h"
namespace mrna {

enum Ctd : int8_t {
  CTD_NA,
  CTD_UNUSED,
  CTD_3_DANGLE,
  CTD_5_DANGLE,
  CTD_MISMATCH,
  CTD_LCOAX_WITH_NEXT,
  CTD_LCOAX_WITH_PREV,
  CTD_RCOAX_WITH_NEXT,
  CTD_RCOAX_WITH_PREV,
  CTD_FCOAX_WITH_NEXT,
  CTD_FCOAX_WITH_PREV,
  CTD_SIZE
};

using Ctds = std::vector<Ctd>;

std::tuple<Primary, Secondary, Ctds> ParsePrimaryCtdString(
    const std::string& prim_str, const std::string& pairs_str);
std::string CtdString(const Secondary& s, const Ctds& ctd);
bool IsCtdString(const std::string& pairs_str);

}  // namespace mrna

#endif  // MODEL_CTD_H_
