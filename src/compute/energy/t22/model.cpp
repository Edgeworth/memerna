// Copyright 2022 Eliot Courtney.
#include "compute/energy/t22/model.h"

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
#include "compute/energy/common/model.h"
#include "compute/energy/common/parse.h"
#include "compute/energy/common/t04like/branch.h"
#include "compute/energy/common/t04like/model_mixin.h"
#include "compute/energy/energy.h"
#include "compute/energy/energy_cfg.h"
#include "compute/energy/structure.h"
#include "model/base.h"
#include "model/ctd.h"
#include "model/primary.h"
#include "util/argparse.h"
#include "util/error.h"
#include "util/string.h"

namespace mrna::erg::t22 {

void Model::LoadFromDir(const std::string& data_dir) {
  T04ModelMixin::LoadFromDir(data_dir);
  Parse4MapFromFile(data_dir + "/penultimate_stacking.data", penultimate_stack);
}

void Model::LoadRandom(std::mt19937& eng) {
  T04ModelMixin::LoadRandom(eng);

  // penultimate_stack is dependent on the direction, so 180 degree rotations
  // don't have to be the same.
  std::uniform_real_distribution<double> energy_dist(RAND_MIN_ENERGY, RAND_MAX_ENERGY);
  RANDOMISE_DATA(penultimate_stack);
}

}  // namespace mrna::erg::t22
