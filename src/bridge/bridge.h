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

const std::map<std::string, ArgParse::option_t> MEMERNA_OPTIONS = {
    {"alg", ArgParse::option_t("which algorithm for memerna").Arg("slow", {"slow", "1", "2", "3"})},
    {"data-path", ArgParse::option_t("path to data tables for memerna").Arg("data/")}
};

// Note that only one energy model can be loaded at a time.
class Memerna : public RnaPackage {
public:
  Memerna(const std::string& data_path, energy_t (* fold_alg_)() = &fold::Fold);
  Memerna(const Memerna&) = delete;
  Memerna& operator=(const Memerna&) = delete;

  virtual energy_t Efn(const folded_rna_t& frna, std::string* desc = nullptr) const;
  virtual folded_rna_t Fold(const rna_t& rna) const;
private:
  energy_t (* fold_alg)();
};

// TODO: change to be able to create all of them
std::unique_ptr<Memerna> MemernaFromArgParse(const ArgParse& argparse);

}
}

#endif //MEMERNA_BRIDGE_H
