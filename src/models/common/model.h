// Copyright 2023 Eliot Courtney.
#ifndef MODELS_COMMON_MODEL_H_
#define MODELS_COMMON_MODEL_H_

#include <memory>
#include <random>
#include <string>

#include "api/energy/energy.h"
#include "api/energy/energy_cfg.h"
#include "model/energy.h"
#include "util/argparse.h"
#include "util/util.h"

namespace mrna::md {

// TODO(1): Get rid of these global macros.

#define CHECK_COND(cond, reason_str)                           \
  do {                                                         \
    if (!(cond)) {                                             \
      if (reason) *reason = "expected " #cond ": " reason_str; \
      return false;                                            \
    }                                                          \
  } while (0)

#define RANDOMISE_DATA(em, d)                                       \
  do {                                                              \
    auto dp = reinterpret_cast<Energy*>(&((em).d));                 \
    /* NOLINTNEXTLINE */                                            \
    for (unsigned int i = 0; i < sizeof((em).d) / sizeof(*dp); ++i) \
      dp[i] = Energy::FromDouble(energy_dist(eng));                 \
  } while (0)

template <typename T>
class ModelMixin {
 public:
  using Ptr = std::shared_ptr<T>;

  static Ptr Create() { return Ptr(new T); }

  static Ptr FromModelPath(const std::string& path) {
    auto em = Create();
    em->LoadFromModelPath(path);
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
    if (args.Has(erg::OPT_SEED)) {
      em = Random(args.Get<uint_fast32_t>(erg::OPT_SEED));
    } else {
      em = FromModelPath(erg::ModelPathFromArgParse(args));
    }
    em->cfg = erg::EnergyCfg::FromArgParse(args);
    return em;
  }

  inline Ptr Clone() const { return std::make_shared<T>(*reinterpret_cast<const T*>(this)); }

  inline Ptr CloneWithPseudofreeEnergy(
      std::vector<Energy> pf_paired, std::vector<Energy> pf_unpaired) {
    auto em = Clone();
    em->LoadPseudofreeEnergy(std::move(pf_paired), std::move(pf_unpaired));
    return em;
  }
};

}  // namespace mrna::md

#endif  // MODELS_COMMON_MODEL_H_
