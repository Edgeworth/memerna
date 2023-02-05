// Copyright 2016 Eliot Courtney.
#include "api/bridge/memerna.h"

#include <memory>
#include <string>
#include <vector>

#include "model/primary.h"
#include "model/structure.h"

namespace mrna::bridge {

erg::EnergyResult Memerna::Efn(const Primary& r, const Secondary& s, std::string* desc) const {
  auto res = ctx_.Efn(r, s, nullptr, desc != nullptr);
  if (desc) {
    for (const auto& struc : res.struc->Description()) {
      *desc += struc;
      *desc += "\n";
    }
  }

  return res;
}

FoldResult Memerna::Fold(const Primary& r) const { return ctx_.Fold(r); }

int Memerna::Suboptimal(subopt::SuboptCallback fn, const Primary& r, Energy delta) const {
  return ctx_.Suboptimal(r, fn, {.delta = delta, .sorted = true});
}

std::vector<subopt::SuboptResult> Memerna::SuboptimalIntoVector(
    const Primary& r, Energy delta) const {
  return ctx_.SuboptimalIntoVector(r, {.delta = delta, .sorted = true});
}

part::PartResult Memerna::Partition(const Primary& r) const { return ctx_.Partition(r); }

Memerna Memerna::FromArgParse(const ArgParse& args) { return Memerna(Ctx::FromArgParse(args)); }

}  // namespace mrna::bridge
