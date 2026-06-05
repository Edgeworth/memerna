// Copyright 2026 Eliot Courtney.
#ifndef BACKENDS_CANON_ENERGY_MODEL_H_
#define BACKENDS_CANON_ENERGY_MODEL_H_

#include <memory>
#include <random>
#include <string>

#include "api/ctx/backend_cfg.h"
#include "api/energy/pseudofree_cfg.h"
#include "backends/base/energy/model.h"
#include "util/util.h"

namespace mrna::md::canon {

class Model : public base::Model {
 public:
  using Ptr = std::shared_ptr<Model>;

  static constexpr auto KIND = BackendKind::CANON;

  static bool IsSupported(const BackendCfg& backend_cfg, const erg::EnergyCfg& cfg,
      const erg::PseudofreeCfg& pf, std::string* reason = nullptr);

  static Ptr Create() { return Ptr(new Model); }

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
    verify(backend == KIND, "expected backend kind: {}, got: {}", KIND, backend);
    return std::visit(
        overloaded{
            [](std::monostate) -> Ptr { fatal("cannot create model without data source"); },
            [&](const std::string&) { return FromModelPath(*cfg.ModelPath(backend)); },
            [&](const RandomModelCfg& random_cfg) { return Random(cfg, random_cfg); },
        },
        cfg.data_src);
  }

  Ptr Clone() const { return std::make_shared<Model>(*this); }

 private:
  // This is private to prevent construction on the stack, since this structure is large.
  Model() = default;
};

}  // namespace mrna::md::canon

#endif  // BACKENDS_CANON_ENERGY_MODEL_H_
