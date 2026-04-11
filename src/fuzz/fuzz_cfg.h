// Copyright 2021 Eliot Courtney.
#ifndef FUZZ_FUZZ_CFG_H_
#define FUZZ_FUZZ_CFG_H_

#include <fmt/ostream.h>

#include <iosfwd>
#include <optional>
#include <string>
#include <vector>

#include "api/ctx/backend_cfg.h"
#include "api/energy/energy_cfg.h"
#include "model/energy.h"
#include "util/argparse.h"

namespace mrna::fuzz {

inline const Opt OPT_FUZZ_BACKENDS = Opt(Opt::ARG)
                                         .LongName("backends")
                                         .ChoiceEnum<BackendKind>()
                                         .Multiple()
                                         .Help("backends to fuzz");

inline const auto OPT_FUZZ_RANDOM_SEEDS =
    mrna::Opt(mrna::Opt::FLAG)
        .LongName("random-seeds")
        .Help("fuzz random energy models with a different seed each time. c.f. seed");

inline const auto OPT_FUZZ_RANDOM_PSEUDOFREE =
    mrna::Opt(mrna::Opt::FLAG).LongName("random-pf").Help("fuzz random pseudofree energies");

inline const auto OPT_FUZZ_RANDOM_PSEUDOFREE_MIN_ENERGY =
    mrna::Opt(mrna::Opt::ARG)
        .LongName("random-pf-min-energy")
        .Default(E(-10.0))
        .Help("minimum energy for generated random pseudofree values");

inline const auto OPT_FUZZ_RANDOM_PSEUDOFREE_MAX_ENERGY =
    mrna::Opt(mrna::Opt::ARG)
        .LongName("random-pf-max-energy")
        .Default(E(10.0))
        .Help("maximum energy for generated random pseudofree values");

// Brute force specific options:
// Allows brute force fuzzing to be given a maximum RNA size
inline const Opt OPT_FUZZ_BRUTE_MAX =
    Opt(Opt::ARG).LongName("brute-max").Default(30).Help("maximum RNA size to run brute force on");

// MFE fuzzing options:
inline const Opt OPT_FUZZ_MFE = Opt(Opt::FLAG).LongName("mfe").Default(0).Help("fuzz mfe");
inline const Opt OPT_FUZZ_MFE_RNASTRUCTURE =
    Opt(Opt::FLAG).LongName("mfe-rnastructure").Default(0).Help("fuzz RNAstructure mfe");
inline const Opt OPT_FUZZ_MFE_TABLE =
    Opt(Opt::FLAG).LongName("mfe-table").Default(0).Help("enable checking mfe dp tables");

// Subopt fuzzing:
inline const Opt OPT_FUZZ_SUBOPT =
    Opt(Opt::FLAG).LongName("subopt").Default(0).Help("fuzz suboptimal folding");
inline const Opt OPT_FUZZ_SUBOPT_RNASTRUCTURE = Opt(Opt::FLAG)
                                                    .LongName("subopt-rnastructure")
                                                    .Default(0)
                                                    .Help("fuzz RNAstructure suboptimal folding");
inline const Opt OPT_FUZZ_SUBOPT_STRUCS =
    Opt(Opt::ARG)
        .LongName("subopt-strucs")
        .Default(10000)
        .Help("maximum number of substructures for subopt fuzz");
inline const Opt OPT_FUZZ_SUBOPT_DELTA = Opt(Opt::ARG)
                                             .LongName("subopt-delta")
                                             .Default(E(0.6))
                                             .Help("max energy delta for subopt delta fuzz");

// Partition fuzzing:
inline const Opt OPT_FUZZ_PFN =
    Opt(Opt::FLAG).LongName("pfn").Default(0).Help("fuzz partition function");
inline const Opt OPT_FUZZ_PFN_RNASTRUCTURE = Opt(Opt::FLAG)
                                                 .LongName("pfn-rnastructure")
                                                 .Default(0)
                                                 .Help("fuzz RNAstructure partition function");
inline const Opt OPT_FUZZ_PFN_SUBOPT =
    Opt(Opt::ARG)
        .LongName("pfn-subopt")
        .Default(0)
        .Help("fuzz partition q value vs suboptimal folding max length");
inline const Opt OPT_FUZZ_PFN_PQ_REL_EP =
    Opt(Opt::ARG)
        .LongName("pfn-pq-rel-ep")
        .Default(EP)
        .Help("partition function Q and P table relative epsilon");
inline const Opt OPT_FUZZ_PFN_PQ_ABS_EP =
    Opt(Opt::ARG)
        .LongName("pfn-pq-abs-ep")
        .Default(EP)
        .Help("partition function Q and P table absolute epsilon");
inline const Opt OPT_FUZZ_PFN_PROB_REL_EP =
    Opt(Opt::ARG)
        .LongName("pfn-prob-rel-ep")
        .Default(EP)
        .Help("partition function probability table relative epsilon");
inline const Opt OPT_FUZZ_PFN_PROB_ABS_EP =
    Opt(Opt::ARG)
        .LongName("pfn-prob-abs-ep")
        .Default(EP)
        .Help("partition function probability table absolute epsilon");

void RegisterOpts(ArgParse* args);

struct RandomPseudofreeCfg {
  RandomPseudofreeCfg(Energy min_energy, Energy max_energy);
  static RandomPseudofreeCfg FromArgParse(const ArgParse& args);

  Energy min_energy;
  Energy max_energy;
};

std::ostream& operator<<(std::ostream& str, const RandomPseudofreeCfg& o);

// Contains all configuration needed to run a fuzzing round.
struct FuzzCfg {
  int brute_max = 22;

  // MFE configuration.
  bool mfe = false;
  bool mfe_rnastructure = false;
  bool mfe_table = false;

  // Subopt folding configuration.
  bool subopt = false;
  bool subopt_rnastructure = false;
  int64_t subopt_strucs = 10000;
  Energy subopt_delta = E(0.6);

  // Partition function configuration.
  bool pfn = false;
  bool pfn_rnastructure = false;
  int pfn_subopt = 0;
  flt pfn_pq_rel_ep = EP;
  flt pfn_pq_abs_ep = EP;
  flt pfn_prob_rel_ep = EP;
  flt pfn_prob_abs_ep = EP;

  // Whether to use a new random-model seed every time.
  bool random_seeds = false;
  bool random_pseudofree = false;

  // Random model range, with an optional fixed seed when random_seeds is false.
  RandomModelCfg random_model_cfg{std::nullopt, E(-10.0), E(10.0)};
  RandomPseudofreeCfg random_pseudofree_cfg{E(-10.0), E(10.0)};

  erg::EnergyCfg energy_cfg{};
  erg::EnergyModelKind energy_model{};
  std::vector<BackendKind> backends{};
  std::vector<Energy> pf_paired{};
  std::vector<Energy> pf_unpaired{};

  std::string data_dir{};
  std::string rnastructure_data_dir{};

  [[nodiscard]] std::string Desc() const;

  static FuzzCfg FromArgParse(const ArgParse& args);
};

}  // namespace mrna::fuzz

template <>
struct fmt::formatter<mrna::fuzz::RandomPseudofreeCfg> : ostream_formatter {};

#endif  // FUZZ_FUZZ_CFG_H_
