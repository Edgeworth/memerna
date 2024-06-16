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

FoldResult Memerna::Fold(const Primary& r) const { return ctx_.Fold(r, {}); }

int Memerna::Subopt(subopt::SuboptCallback fn, const Primary& r, Energy delta) const {
  return ctx_.Subopt(r, fn, {.delta = delta, .sorted = true});
}

std::vector<subopt::SuboptResult> Memerna::SuboptIntoVector(const Primary& r, Energy delta) const {
  return ctx_.SuboptIntoVector(r, {.delta = delta, .sorted = true});
}

pfn::PfnResult Memerna::Pfn(const Primary& r) const { return ctx_.Pfn(r); }

Memerna Memerna::FromArgParse(const ArgParse& args) { return Memerna(Ctx::FromArgParse(args)); }

}  // namespace mrna::bridge
