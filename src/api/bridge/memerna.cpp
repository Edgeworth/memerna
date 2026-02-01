// Copyright 2016 Eliot Courtney.
#include "api/bridge/memerna.h"

#include <memory>
#include <string>
#include <vector>

#include "api/trace/trace_cfg.h"
#include "model/primary.h"
#include "model/structure.h"

namespace mrna::bridge {

erg::EnergyResult Memerna::Efn(const Primary& r, const Secondary& s, std::string* desc) const {
  // TODO(2): Support pseudofree energy in the bridge API.
  auto res = ctx_.Efn(r, s, erg::EnergyCfg{}, erg::PseudofreeCfg{}, /*given_ctd=*/nullptr,
      /*build_structure=*/desc != nullptr);
  if (desc) {
    for (const auto& struc : res.struc->Description()) {
      *desc += struc;
      *desc += "\n";
    }
  }

  return res;
}

FoldResult Memerna::Fold(const Primary& r) const {
  return ctx_.Fold(r, erg::EnergyCfg{}, erg::PseudofreeCfg{}, trace::TraceCfg{});
}

int Memerna::Subopt(subopt::SuboptCallback fn, const Primary& r, Energy delta) const {
  return ctx_.Subopt(
      r, erg::EnergyCfg{}, erg::PseudofreeCfg{}, fn, {.delta = delta, .sorted = true});
}

std::vector<subopt::SuboptResult> Memerna::SuboptIntoVector(const Primary& r, Energy delta) const {
  return ctx_.SuboptIntoVector(
      r, erg::EnergyCfg{}, erg::PseudofreeCfg{}, {.delta = delta, .sorted = true});
}

pfn::PfnResult Memerna::Pfn(const Primary& r) const {
  return ctx_.Pfn(r, erg::EnergyCfg{}, erg::PseudofreeCfg{});
}

Memerna Memerna::FromArgParse(const ArgParse& args) { return Memerna(Ctx::FromArgParse(args)); }

}  // namespace mrna::bridge
