#ifndef MEMERNA_BRIDGE_H
#define MEMERNA_BRIDGE_H

#include "argparse.h"
#include "common.h"
#include "fold/fold.h"

#include "nn_unpaired_model.hpp"
#include "RNAstructure/rna_library.h"

namespace memerna {
namespace bridge {

class RnaPackage {
public:

  virtual ~RnaPackage() = default;

  virtual energy_t Efn(const folded_rna_t& frna, std::string* desc = nullptr) const = 0;

  virtual folded_rna_t Fold(const rna_t& rna) const = 0;
};

class Rnark : public RnaPackage {
public:
  Rnark(const std::string& data_path);
  Rnark(const Rnark&) = delete;
  Rnark& operator=(const Rnark&) = delete;

  virtual energy_t Efn(const folded_rna_t& frna, std::string* desc = nullptr) const;
  virtual folded_rna_t Fold(const rna_t& rna) const;

private:
  librnary::NNUnpairedModel model;
};

class Rnastructure : public RnaPackage {
public:
  Rnastructure(const std::string& data_path, bool use_lyngso_);
  Rnastructure(const Rnastructure&) = delete;
  Rnastructure& operator=(Rnastructure&) = delete;

  virtual energy_t Efn(const folded_rna_t& frna, std::string* desc = nullptr) const;
  virtual folded_rna_t Fold(const rna_t& rna) const;

private:
  std::unique_ptr<datatable> data;
  bool use_lyngso;
};

// Note that only one energy model can be loaded at a time.
class Memerna : public RnaPackage {
public:
  Memerna(const std::string& data_path, const fold::fold_fn_t* fold_fn_ = &fold::Fold0);
  Memerna(const Memerna&) = delete;
  Memerna& operator=(const Memerna&) = delete;

  virtual energy_t Efn(const folded_rna_t& frna, std::string* desc = nullptr) const;
  virtual folded_rna_t Fold(const rna_t& rna) const;
private:
  fold::fold_fn_t* fold_fn;
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
