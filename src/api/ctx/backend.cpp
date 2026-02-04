// Copyright 2024 Eliot Courtney.
#include "api/ctx/backend.h"

#include "util/error.h"

namespace mrna {

BackendModelPtr BackendFromArgParse(const ArgParse& args) {
  return BackendFromBackendCfg(BackendCfg::FromArgParse(args));
}

BackendModelPtr BackendFromBackendCfg(const BackendCfg& cfg) {
  BackendCfg resolved = cfg;
  if (resolved.backend == BackendKind::AUTO) resolved.backend = BackendKind::BASEOPT;
  switch (resolved.backend) {
  case BackendKind::AUTO: break;
  case BackendKind::BASE: return md::base::Model::FromBackendCfg(resolved);
  case BackendKind::BASEOPT: return md::base::opt::Model::FromBackendCfg(resolved);
  case BackendKind::STACK: return md::stack::Model::FromBackendCfg(resolved);
  }
  unreachable();
}

BackendModelPtr Random(BackendKind kind, uint_fast32_t seed) {
  if (kind == BackendKind::AUTO) kind = BackendKind::BASEOPT;
  switch (kind) {
  case BackendKind::AUTO: break;
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
