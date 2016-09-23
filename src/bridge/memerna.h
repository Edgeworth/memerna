#ifndef MEMERNA_MEMERNA_H
#define MEMERNA_MEMERNA_H

#include "common.h"
#include "bridge/bridge.h"

namespace memerna {
namespace bridge {

// Note that only one energy model can be loaded at a time.
class Memerna : public RnaPackage {
public:
  Memerna(energy::EnergyModelPtr em_, fold::context_options_t options_) : em(em_), options(options_) {}

  Memerna(const Memerna&) = delete;
  Memerna& operator=(const Memerna&) = delete;

  virtual energy_t Efn(const secondary_t& secondary, std::string* desc = nullptr) const;
  virtual computed_t Fold(const primary_t& r) const;
  virtual std::vector<computed_t> Suboptimal(const primary_t& r, energy_t energy_delta) const;
private:
  const energy::EnergyModelPtr em;
  const fold::context_options_t options;
};


}
}

#endif  // MEMERNA_MEMERNA_H
