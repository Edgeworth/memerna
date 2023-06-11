// Copyright 2022 Eliot Courtney.
#include "api/energy/energy.h"

#include <utility>

#include "api/energy/energy_cfg.h"
#include "model/structure.h"

namespace mrna::erg {

void RegisterOptsEnergyModel(ArgParse* args) {
  RegisterOptsEnergyCfg(args);
  args->RegisterOpt(OPT_MEMERNA_DATA);
  args->RegisterOpt(OPT_ENERGY_MODEL);
  args->RegisterOpt(OPT_SEED);
  args->RegisterOpt(OPT_PAIRED_PSEUDOFREE);
  args->RegisterOpt(OPT_UNPAIRED_PSEUDOFREE);
}

std::string ModelPathFromArgParse(const ArgParse& args) {
  return ModelPathFromArgParse(args, args.Get(OPT_ENERGY_MODEL));
}

std::string ModelPathFromArgParse(const ArgParse& args, const std::string& model) {
  return ModelPath(args.Get(OPT_MEMERNA_DATA), model);
}

std::string ModelPath(const std::string& data_dir, const std::string& model) {
  return data_dir + "/model/" + model;
}

EnergyResult::EnergyResult(Energy energy, Ctds ctd, std::unique_ptr<Structure> struc)
    : energy(energy), ctd(std::move(ctd)), struc(std::move(struc)) {}

EnergyResult::~EnergyResult() = default;

}  // namespace mrna::erg
