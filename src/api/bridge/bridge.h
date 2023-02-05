// Copyright 2016 Eliot Courtney.
#ifndef BRIDGE_BRIDGE_H_
#define BRIDGE_BRIDGE_H_

#include <memory>
#include <string>
#include <vector>

#include "api/part.h"
#include "api/subopt/subopt.h"
#include "api/ctx/ctx.h"
#include "model/energy.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "util/argparse.h"

namespace mrna::bridge {

inline const Opt OPT_USE_RNASTRUCTURE =
    Opt(Opt::FLAG).LongName("rnastructure").ShortName("r").Default(false);
inline const Opt OPT_RNASTRUCTURE_DATA =
    Opt(Opt::ARG).LongName("rnastructure-data").ShortName("rd").Help("data path for RNAstructure");
inline const Opt OPT_USE_MEMERNA = Opt(Opt::FLAG).LongName("memerna").ShortName("m").Default(false);

void RegisterOpts(ArgParse* args);

class RnaPackage {
 public:
  RnaPackage() = default;
  virtual ~RnaPackage() = default;

  RnaPackage(RnaPackage&& o) = default;
  RnaPackage& operator=(RnaPackage&&) = default;

  RnaPackage(const RnaPackage&) = delete;
  RnaPackage& operator=(const RnaPackage&) = delete;

  virtual erg::EnergyResult Efn(
      const Primary& r, const Secondary& s, std::string* desc = nullptr) const = 0;
  [[nodiscard]] virtual ctx::FoldResult Fold(const Primary& r) const = 0;
  [[nodiscard]] virtual int Suboptimal(
      subopt::SuboptCallback fn, const Primary& r, Energy delta) const = 0;
  [[nodiscard]] virtual std::vector<subopt::SuboptResult> SuboptimalIntoVector(
      const Primary& r, Energy delta) const = 0;
  [[nodiscard]] virtual part::PartResult Partition(const Primary& r) const = 0;

  static std::unique_ptr<RnaPackage> FromArgParse(const ArgParse& args);
};

}  // namespace mrna::bridge

#endif  // BRIDGE_BRIDGE_H_
