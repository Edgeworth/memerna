// Copyright 2016 Eliot Courtney.
#include "compute/energy/t04/model.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "compute/energy/common/branch.h"
#include "compute/energy/common/t04like/branch.h"
#include "compute/energy/energy.h"
#include "compute/energy/energy_cfg.h"
#include "compute/energy/structure.h"
#include "model/base.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "util/argparse.h"
#include "util/error.h"
#include "util/string.h"

namespace mrna::erg::t04 {}  // namespace mrna::erg::t04
