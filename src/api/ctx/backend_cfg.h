// Copyright 2024 Eliot Courtney.
#ifndef API_CTX_BACKEND_CFG_H_
#define API_CTX_BACKEND_CFG_H_

#include <fmt/ostream.h>

#include <cstdint>
#include <iosfwd>
#include <optional>
#include <string>
#include <variant>

#include "api/energy/energy_cfg.h"
#include "model/energy.h"
#include "util/argparse.h"
#include "util/string.h"

namespace mrna {

void RegisterOptsBackendCfg(ArgParse* args);

MAKE_ENUM(BackendKind, BASE, BASEOPT, STACK);

struct RandomModelCfg {
  RandomModelCfg(std::optional<uint_fast32_t> seed, Energy min_energy, Energy max_energy);
  static RandomModelCfg FromArgParse(const ArgParse& args);

  std::optional<uint_fast32_t> seed;
  Energy min_energy;
  Energy max_energy;
};

std::ostream& operator<<(std::ostream& str, const RandomModelCfg& o);

struct BackendCfg {
  erg::EnergyModelKind energy_model = erg::EnergyModelKind::T04;
  int precision = MRNA_ENERGY_PRECISION;
  std::variant<std::monostate, std::string, RandomModelCfg> data_src = std::monostate{};

  [[nodiscard]]
  static BackendCfg FromArgParse(const ArgParse& args);

  // Returns the model path if data_src is a data directory, nullopt otherwise.
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
                                            .Choice({Conv(MRNA_ENERGY_PRECISION)})
                                            .Default(MRNA_ENERGY_PRECISION)
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

inline const Opt OPT_RAND_MIN_ENERGY = Opt(Opt::ARG)
                                           .LongName("random-min-energy")
                                           .Default(E(-10.0))
                                           .Help("minimum energy for generated random models");

inline const Opt OPT_RAND_MAX_ENERGY = Opt(Opt::ARG)
                                           .LongName("random-max-energy")
                                           .Default(E(10.0))
                                           .Help("maximum energy for generated random models");

}  // namespace mrna

template <>
struct fmt::formatter<mrna::RandomModelCfg> : ostream_formatter {};

#endif  // API_CTX_BACKEND_CFG_H_
