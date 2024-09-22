#include "backends/common/base/model_base.h"

namespace mrna::md::base {

// static
void ModelBase::VerifyForEfn(const Primary& r, const Secondary& s, const Ctds* given_ctd) {
  verify(r.size() == s.size(), "sequence and secondary structure must be the same length");
  if (given_ctd)
    verify(given_ctd->size() == r.size(), "given CTDs must be the same length as the seq");
}

}  // namespace mrna::md::base
