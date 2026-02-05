// Copyright 2024 Eliot Courtney.
#include "api/ctx/backend.h"

#include "util/error.h"

namespace mrna {

BackendModelPtr BackendFromBackendCfg(BackendKind backend, const BackendCfg& cfg) {
  switch (backend) {
  case BackendKind::BASE: return md::base::Model::FromBackendCfg(backend, cfg);
  case BackendKind::BASEOPT: return md::base::opt::Model::FromBackendCfg(backend, cfg);
  case BackendKind::STACK: return md::stack::Model::FromBackendCfg(backend, cfg);
  }
  unreachable();
}

BackendModelPtr Random(BackendKind kind, uint_fast32_t seed) {
  switch (kind) {
  case BackendKind::BASE: return md::base::Model::Random(seed);
  case BackendKind::BASEOPT: return md::base::opt::Model::Random(seed);
  case BackendKind::STACK: return md::stack::Model::Random(seed);
  }
  unreachable();
}

BackendBoltzModelPtr Boltz(const BackendModelPtr& m) {
  auto vis = overloaded{
      [](const md::base::Model::Ptr& m) -> BackendBoltzModelPtr {
        return md::base::BoltzModel::Create(m);
      },
      [](const md::base::opt::Model::Ptr& m) -> BackendBoltzModelPtr {
        return md::base::opt::BoltzModel::Create(m);
      },
      [](const md::stack::Model::Ptr& m) -> BackendBoltzModelPtr {
        return md::stack::BoltzModel::Create(m);
      },
  };
  return std::visit(vis, m);
}

BackendModelPtr CloneBackend(const BackendModelPtr& m) {
  auto vis = overloaded{[](const auto& m) { return BackendModelPtr(m->Clone()); }};
  return std::visit(vis, m);
}

}  // namespace mrna
