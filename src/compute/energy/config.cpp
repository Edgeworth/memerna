// Copyright 2022 Eliot Courtney.
#include "compute/energy/config.h"

namespace mrna::energy {

void RegisterOpts(ArgParse* args) {
  args->RegisterOpt(OPT_SEED);
  args->RegisterOpt(OPT_MEMERNA_DATA);
}

}  // namespace mrna::energy
