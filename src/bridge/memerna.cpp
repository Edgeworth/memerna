// Copyright 2016 Eliot Courtney.
#include "bridge/memerna.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "compute/energy/structure.h"
#include "model/primary.h"

namespace mrna::bridge {

energy::EnergyResult Memerna::Efn(const Primary& r, const Secondary& s, std::string* desc) const {
  auto res = ctx_.Efn(r, s, nullptr, desc != nullptr);
  if (desc) {
    for (const auto& s : res.struc->Description()) {
      *desc += s;
      *desc += "\n";
    }
  }

  return res;
}

ctx::FoldResult Memerna::Fold(const Primary& r) const { return ctx_.Fold(r); }

int Memerna::Suboptimal(subopt::SuboptCallback fn, const Primary& r, Energy delta) const {
  return ctx_.Suboptimal(r, fn, {.delta = delta, .sorted = true});
}

std::vector<subopt::SuboptResult> Memerna::SuboptimalIntoVector(
    const Primary& r, Energy delta) const {
  return ctx_.SuboptimalIntoVector(r, {.delta = delta, .sorted = true});
}

part::PartResult Memerna::Partition(const Primary& r) const { return ctx_.Partition(r); }

Memerna Memerna::FromArgParse(const ArgParse& args) {
  return Memerna(ctx::Ctx::FromArgParse(args));
}

}  // namespace mrna::bridge
