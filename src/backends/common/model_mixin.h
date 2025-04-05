// Copyright 2023 Eliot Courtney.
#ifndef BACKENDS_COMMON_MODEL_MIXIN_H_
#define BACKENDS_COMMON_MODEL_MIXIN_H_

#include <memory>
#include <random>
#include <string>

#include "api/ctx/backend_cfg.h"

namespace mrna::md {

// TODO(1): Get rid of these global macros.

#define CHECK_COND(cond, reason_str)                           \
  do {                                                         \
    if (!(cond)) {                                             \
      if (reason) *reason = "expected " #cond ": " reason_str; \
      return false;                                            \
    }                                                          \
  } while (0)

#define RANDOMISE_DATA(m, d)                                       \
  do {                                                             \
    auto dp = reinterpret_cast<Energy*>(&((m).d));                 \
    /* NOLINTNEXTLINE */                                           \
    for (unsigned int i = 0; i < sizeof((m).d) / sizeof(*dp); ++i) \
      dp[i] = Energy::FromFlt(energy_dist(eng));                   \
  } while (0)

template <typename T>
class ModelMixin {
 public:
  using Ptr = std::shared_ptr<T>;

  static Ptr Create() { return Ptr(new T); }

  static Ptr FromModelPath(const std::string& path) {
    auto m = Create();
    m->LoadFromModelPath(path);
    std::string reason;
    verify(m->IsValid(&reason), "invalid energy model: {}", reason);
    return m;
  }

  static Ptr Random(uint_fast32_t seed) {
    auto m = Create();
    std::mt19937 eng(seed);
    m->LoadRandom(eng);
    std::string reason;
    verify(m->IsValid(&reason), "invalid energy model: {}", reason);
    return m;
  }

  static Ptr FromBackendCfg(const BackendCfg& cfg) {
    verify(cfg.backend == T::KIND, "expected backend kind: {}, got: {}", T::KIND, cfg.backend);
    Ptr m;
    if (cfg.seed.has_value()) {
      m = Random(*cfg.seed);
    } else {
      m = FromModelPath(cfg.BackendDataPath());
    }
    if (!cfg.pf_paired.empty() || !cfg.pf_unpaired.empty()) {
      m->LoadPseudofreeEnergy(cfg.pf_paired, cfg.pf_unpaired);
    }
    m->SetEnergyCfg(cfg.energy_cfg);
    return m;
  }

  Ptr Clone() const { return std::make_shared<T>(*reinterpret_cast<const T*>(this)); }
};

}  // namespace mrna::md

#endif  // BACKENDS_COMMON_MODEL_MIXIN_H_
