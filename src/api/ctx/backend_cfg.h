// Copyright 2024 Eliot Courtney.
#ifndef API_CTX_BACKEND_CFG_H_
#define API_CTX_BACKEND_CFG_H_

#include <cstdint>
#include <vector>

#include "api/energy/energy_cfg.h"
#include "model/energy.h"
#include "util/argparse.h"
#include "util/string.h"

namespace mrna {

void RegisterOptsBackendCfg(ArgParse* args);

MAKE_ENUM(BackendKind, BASE, BASEOPT, STACK);

struct BackendCfg {
  erg::EnergyModelKind energy_model = erg::EnergyModelKind::T04;
  int precision = ENERGY_PRECISION;
  BackendKind backend = BackendKind::BASEOPT;
  std::string data_dir{};
  std::optional<uint_fast32_t> seed = std::nullopt;
  erg::EnergyCfg energy_cfg{};
  std::vector<Energy> pf_paired{};
  std::vector<Energy> pf_unpaired{};

  [[nodiscard]]
  static BackendCfg FromArgParse(const ArgParse& args);

  [[nodiscard]]
  std::string BackendDataPath() const;
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

inline const Opt OPT_BACKEND = Opt(Opt::ARG)
                                   .LongName("backend")
                                   .ShortName("b")
                                   .ChoiceEnum<BackendKind>()
                                   .Default(BackendKind::BASEOPT)
                                   .Help("backend to use");

inline const Opt OPT_MEMERNA_DATA = Opt(Opt::ARG)
                                        .LongName("memerna-data")
                                        .ShortName("md")
                                        .Default("./data/")
                                        .Help("data path for memerna data");

inline const Opt OPT_SEED =
    Opt(Opt::ARG).LongName("seed").Help("seed for random energy model for memerna");

inline const Opt OPT_PAIRED_PSEUDOFREE =
    Opt(Opt::ARG)
        .LongName("pf-paired")
        .Multiple()
        .Help("comma separated energies for paired pseudofree energy");
inline const Opt OPT_UNPAIRED_PSEUDOFREE =
    Opt(Opt::ARG)
        .LongName("pf-unpaired")
        .Multiple()
        .Help("comma separated energies for unpaired pseudofree energy");

}  // namespace mrna

#endif  // API_CTX_BACKEND_CFG_H_
