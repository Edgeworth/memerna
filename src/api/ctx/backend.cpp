// Copyright 2024 Eliot Courtney.
#include "api/ctx/backend.h"

#include <variant>

#include "api/ctx/backend_cfg.h"
#include "backends/base/energy/model.h"
#include "backends/canon/energy/model.h"
#include "backends/stack/energy/model.h"
#include "util/error.h"
#include "util/util.h"

namespace mrna {

BackendModelPtr BackendFromBackendCfg(BackendKind backend, const BackendCfg& cfg) {
  switch (backend) {
  case BackendKind::BASE: return md::base::Model::FromBackendCfg(backend, cfg);
  case BackendKind::STACK: return md::stack::Model::FromBackendCfg(backend, cfg);
  case BackendKind::CANON: return md::canon::Model::FromBackendCfg(backend, cfg);
  }
  unreachable();
}

BackendBoltzModelPtr Boltz(const BackendModelPtr& m) {
  return std::visit([](const auto& m) -> BackendBoltzModelPtr { return Boltz(m); }, m);
}

BackendModelPtr CloneBackend(const BackendModelPtr& m) {
  auto vis = overloaded{[](const auto& m) { return BackendModelPtr(m->Clone()); }};
  return std::visit(vis, m);
}

}  // namespace mrna
