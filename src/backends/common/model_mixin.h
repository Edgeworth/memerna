// Copyright 2023 Eliot Courtney.
#ifndef BACKENDS_COMMON_MODEL_MIXIN_H_
#define BACKENDS_COMMON_MODEL_MIXIN_H_

#include <algorithm>  // IWYU pragma: keep - used by RANDOMISE_DATA.
#include <cstddef>  // IWYU pragma: keep - used by RANDOMISE_DATA.
#include <memory>
#include <random>
#include <string>
#include <variant>

#include "api/ctx/backend_cfg.h"
#include "util/error.h"
#include "util/util.h"

namespace mrna::md {

#define CHECK_COND(cond, reason_str)                           \
  do {                                                         \
    if (!(cond)) {                                             \
      if (reason) *reason = "expected " #cond ": " reason_str; \
      return false;                                            \
    }                                                          \
  } while (0)

// Draws table `d` from a random sub-range of `energy_dist`. Pruning bounds
// that add up the minima of several tables, like those in PrecompBase, are
// only tested when the tables have different minima.
#define RANDOMISE_DATA(m, d)                                                 \
  do {                                                                       \
    auto* data = Decay((m).d);                                               \
    const auto [lo, hi] = std::minmax({energy_dist(eng), energy_dist(eng)}); \
    std::uniform_real_distribution<double> table_dist(lo, hi);               \
    for (std::size_t i = 0; i < sizeof((m).d) / sizeof(*data); ++i)          \
      data[i] = Energy::FromFlt(table_dist(eng));                            \
  } while (0)

template <typename T>
class ModelMixin {
 public:
  using Ptr = std::shared_ptr<T>;

  bool operator==(const ModelMixin&) const = default;

  static Ptr Create() { return Ptr(new T); }

  static Ptr FromModelPath(const std::string& path) {
    auto m = Create();
    m->LoadFromModelPath(path);
    std::string reason;
    verify(m->IsValid(&reason), "invalid energy model: {}", reason);
    return m;
  }

  static Ptr Random(const BackendCfg& cfg, const RandomModelCfg& random_cfg) {
    auto m = Create();
    verify(random_cfg.seed.has_value(), "random model config missing seed");
    std::mt19937 eng(*random_cfg.seed);
    m->LoadRandom(cfg, random_cfg, eng);
    std::string reason;
    verify(m->IsValid(&reason), "invalid energy model: {}", reason);
    return m;
  }

  static Ptr FromBackendCfg(BackendKind backend, const BackendCfg& cfg) {
    verify(backend == T::KIND, "expected backend kind: {}, got: {}", T::KIND, backend);
    return std::visit(
        overloaded{
            [](std::monostate) -> Ptr { fatal("cannot create model without data source"); },
            [&](const std::string&) { return FromModelPath(*cfg.ModelPath(backend)); },
            [&](const RandomModelCfg& random_cfg) { return Random(cfg, random_cfg); },
        },
        cfg.data_src);
  }

  Ptr Clone() const { return std::make_shared<T>(*static_cast<const T*>(this)); }
};

}  // namespace mrna::md

#endif  // BACKENDS_COMMON_MODEL_MIXIN_H_
