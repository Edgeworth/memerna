#ifndef MEMERNA_BRIDGE_H
#define MEMERNA_BRIDGE_H

#include "common.h"

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
  Rnastructure(const std::string& data_path, bool _use_lyngso);
  Rnastructure(const Rnastructure&) = delete;
  Rnastructure& operator=(Rnastructure&) = delete;

  virtual energy_t Efn(const folded_rna_t& frna, std::string* desc = nullptr) const;
  virtual folded_rna_t Fold(const rna_t& rna) const;
private:
  std::unique_ptr<datatable> data;
  bool use_lyngso;
};

class Memerna : public RnaPackage {
public:
  Memerna(const std::string& data_path);
  Memerna(const Memerna&) = delete;
  Memerna& operator=(const Memerna&) = delete;

  virtual energy_t Efn(const folded_rna_t& frna, std::string* desc = nullptr) const;
  virtual folded_rna_t Fold(const rna_t& rna) const;
};

}
}

#endif //MEMERNA_BRIDGE_H
