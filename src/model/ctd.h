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

// Secondary structure with MFE information and CTDs.
struct Computed {
  Computed() = default;
  Computed(const Computed&) = default;
  Computed(Computed&&) = default;
  Computed& operator=(Computed&&) = default;

  explicit Computed(Primary r_) : s(std::move(r_)), base_ctds(s.r.size(), CTD_NA), energy(MAX_E) {}

  explicit Computed(Secondary s_) : s(std::move(s_)), base_ctds(s.r.size(), CTD_NA), energy(MAX_E) {
    verify(s.r.size() == s.p.size() && s.r.size() == base_ctds.size(), "bug");
  }

  Computed(Secondary s_, std::vector<Ctd> base_ctds_, Energy energy_)
      : s(std::move(s_)), base_ctds(std::move(base_ctds_)), energy(energy_) {
    verify(s.r.size() == s.p.size() && s.r.size() == base_ctds.size(), "bug");
  }

  bool operator==(const Computed& o) const {
    return s == o.s && base_ctds == o.base_ctds && energy == o.energy;
  }
  bool operator!=(const Computed& o) const { return !(*this == o); }

  Secondary s;
  std::vector<Ctd> base_ctds;
  Energy energy;
};

struct ComputedEnergyCmp {
  bool operator()(const Computed& a, const Computed& b) const { return a.energy < b.energy; }
};

}  // namespace mrna

#endif  // MODEL_CTD_H_
