// Copyright 2016 Eliot Courtney.
#ifndef BRIDGE_BRIDGE_H_
#define BRIDGE_BRIDGE_H_

#include <memory>
#include <string>
#include <vector>

#include "compute/energy/energy.h"
#include "compute/partition/partition.h"
#include "compute/subopt/subopt.h"
#include "ctx/ctx.h"
#include "model/model.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/argparse.h"

namespace mrna::bridge {

inline const Opt OPT_USE_RNASTRUCTURE = Opt().LongName("rnastructure").ShortName("r");
inline const Opt OPT_RNASTRUCTURE_DATA =
    Opt().LongName("rnastructure-data").ShortName("rd").Arg().Help("data path for RNAstructure");
inline const Opt OPT_USE_MEMERNA = Opt().LongName("memerna").ShortName("m");

void RegisterOpts(ArgParse* args);

class RnaPackage {
 public:
  virtual ~RnaPackage() = default;

  virtual energy::EnergyResult Efn(
      const Primary& r, const Secondary& s, std::string* desc = nullptr) const = 0;
  virtual ctx::FoldResult Fold(const Primary& r) const = 0;
  virtual int Suboptimal(subopt::SuboptCallback fn, const Primary& r, Energy delta) const = 0;
  virtual std::vector<subopt::SuboptResult> SuboptimalIntoVector(
      const Primary& r, Energy delta) const = 0;
  virtual part::PartResult Partition(const Primary& r) const = 0;

  static std::unique_ptr<RnaPackage> FromArgParse(const ArgParse& args);
};

}  // namespace mrna::bridge

#endif  // BRIDGE_BRIDGE_H_
