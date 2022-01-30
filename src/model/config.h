// Copyright 2022 Eliot Courtney.
#ifndef MODEL_CONFIG_H_
#define MODEL_CONFIG_H_

#include "util/argparse.h"

namespace mrna {

// void RegisterOpts(ArgParse* args);

struct ModelCfg {
  bool mfe;
  bool subopt;
  bool part;

  static ModelCfg FromArgParse(const ArgParse& args);
};

}  // namespace mrna

#endif  // MODEL_CONFIG_H_
