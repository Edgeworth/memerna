// Copyright 2023 Eliot Courtney.
#ifndef COMPUTE_ENERGY_COMMON_MODEL_H_
#define COMPUTE_ENERGY_COMMON_MODEL_H_

#include <memory>
#include <random>
#include <string>

#include "compute/energy/energy_cfg.h"
#include "model/energy.h"
#include "util/argparse.h"
#include "util/util.h"

namespace mrna::erg {

// TODO(1): Get rid of these global macros.

#define CHECK_COND(cond, reason_str)                           \
  do {                                                         \
    if (!(cond)) {                                             \
      if (reason) *reason = "expected " #cond ": " reason_str; \
      return false;                                            \
    }                                                          \
  } while (0)

#define RANDOMISE_DATA(d)                                      \
  do {                                                         \
    auto dp = reinterpret_cast<Energy*>(&(d));                 \
    /* NOLINTNEXTLINE */                                       \
    for (unsigned int i = 0; i < sizeof(d) / sizeof(*dp); ++i) \
      dp[i] = Energy::FromDouble(energy_dist(eng));            \
  } while (0)

template <typename T>
class ModelMixin {
 public:
  using Ptr = std::shared_ptr<T>;

  static Ptr Create() { return Ptr(new T); }

  static Ptr FromDir(const std::string& data_dir) {
    auto em = Create();
    em->LoadFromDir(data_dir);
    std::string reason;
    verify(em->IsValid(&reason), "invalid energy model: {}", reason);
    return em;
  }

  static Ptr Random(uint_fast32_t seed) {
    auto em = Create();
    std::mt19937 eng(seed);
    em->LoadRandom(eng);
    std::string reason;
    verify(em->IsValid(&reason), "invalid energy model: {}", reason);
    return em;
  }

  static Ptr FromArgParse(const ArgParse& args) {
    Ptr em;
    if (args.Has(OPT_SEED)) {
      em = Random(args.Get<uint_fast32_t>(OPT_SEED));
    } else {
      em = FromDir(ModelPathFromArgParse(args));
    }
    em->cfg = EnergyCfg::FromArgParse(args);
    return em;
  }

  inline Ptr Clone() const { return std::make_shared<T>(*reinterpret_cast<const T*>(this)); }
};

}  // namespace mrna::erg

#endif  // COMPUTE_ENERGY_COMMON_MODEL_H_
