// Copyright 2024 Eliot Courtney.
#ifndef API_CTX_BACKEND_CFG_H_
#define API_CTX_BACKEND_CFG_H_

#include <cstdint>
#include <optional>
#include <string>
#include <variant>

#include "api/energy/energy_cfg.h"
#include "util/argparse.h"
#include "util/string.h"

namespace mrna {

void RegisterOptsBackendCfg(ArgParse* args);

MAKE_ENUM(BackendKind, BASE, BASEOPT, STACK);

struct BackendCfg {
  erg::EnergyModelKind energy_model = erg::EnergyModelKind::T04;
  int precision = ENERGY_PRECISION;
  std::variant<std::monostate, std::string, uint_fast32_t> data_src = std::monostate{};

  [[nodiscard]]
  static BackendCfg FromArgParse(const ArgParse& args);

  // Returns the model path if data_src is a data directory, nullopt if it's a seed.
  [[nodiscard]]
  std::optional<std::string> ModelPath(BackendKind backend) const;
};

inline const Opt OPT_ENERGY_MODEL = Opt(Opt::ARG)
                                        .LongName("energy-model")
                                        .ShortName("em")
                                        .ChoiceEnum<erg::EnergyModelKind>()
                                        .Default(erg::EnergyModelKind::T04)
                                        .Help("energy model to use");

inline const Opt OPT_ENERGY_PRECISION = Opt(Opt::ARG)
                                            .LongName("energy-precision")
                                            .Choice({Conv(ENERGY_PRECISION)})
                                            .Default(ENERGY_PRECISION)
                                            .Help("energy precision to use");

inline const Opt OPT_BACKEND =
    Opt(Opt::ARG).LongName("backend").ShortName("b").ChoiceEnum<BackendKind>().Help(
        "backend to use (default: auto-select best)");

inline const Opt OPT_MEMERNA_DATA = Opt(Opt::ARG)
                                        .LongName("memerna-data")
                                        .ShortName("md")
                                        .Default("./data/")
                                        .Help("data path for memerna data");

inline const Opt OPT_SEED =
    Opt(Opt::ARG).LongName("seed").Help("seed for random energy model for memerna");

}  // namespace mrna

#endif  // API_CTX_BACKEND_CFG_H_
