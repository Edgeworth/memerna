#ifndef MEMERNA_BRIDGE_H
#define MEMERNA_BRIDGE_H

#include "argparse.h"
#include "common.h"
#include "fold/context.h"

#include "nn_unpaired_model.hpp"
#include "RNAstructure/rna_library.h"

namespace memerna {
namespace bridge {

class RnaPackage {
public:
  virtual ~RnaPackage() = default;

  virtual energy_t Efn(const secondary_t& secondary, std::string* desc = nullptr) const = 0;
  virtual computed_t Fold(const primary_t& r) const = 0;
};

class Rnark : public RnaPackage {
public:
  Rnark(const std::string& data_path);
  Rnark(const Rnark&) = delete;
  Rnark& operator=(const Rnark&) = delete;

  virtual energy_t Efn(const secondary_t& secondary, std::string* desc = nullptr) const;
  virtual computed_t Fold(const primary_t& r) const;

private:
  librnary::NNUnpairedModel model;
};

class Rnastructure : public RnaPackage {
public:
  Rnastructure(const std::string& data_path, bool use_lyngso_);
  Rnastructure(const Rnastructure&) = delete;
  Rnastructure& operator=(Rnastructure&) = delete;

  virtual energy_t Efn(const secondary_t& secondary, std::string* desc = nullptr) const;
  virtual computed_t Fold(const primary_t& r) const;

  computed_t FoldAndDpTable(const primary_t& r, dp_state_t* dp_state) const;

private:
  const std::unique_ptr<datatable> data;
  const bool use_lyngso;
};

// Note that only one energy model can be loaded at a time.
class Memerna : public RnaPackage {
public:
  Memerna(energy::EnergyModelPtr em_, fold::context_options_t options_) : em(em_), options(options_) {}

  Memerna(const Memerna&) = delete;
  Memerna& operator=(const Memerna&) = delete;

  virtual energy_t Efn(const secondary_t& secondary, std::string* desc = nullptr) const;
  virtual computed_t Fold(const primary_t& r) const;
private:
  const energy::EnergyModelPtr em;
  const fold::context_options_t options;
};

const std::map<std::string, ArgParse::option_t> BRIDGE_OPTIONS = {
    {"r", {"rnastructure"}},
    {"m", {"rnark"}},
    {"k", {"memerna"}}
};

std::unique_ptr<RnaPackage> RnaPackageFromArgParse(const ArgParse& argparse);

}
}

#endif //MEMERNA_BRIDGE_H
