#ifndef MEMERNA_FOLD_H
#define MEMERNA_FOLD_H

#include "common.h"
#include "energy.h"

namespace memerna {
namespace fold {

// This function is not re-entrant?
energy_t Fold(std::unique_ptr<structure::Structure>* s = nullptr);

inline energy_t Fold(const rna_t& rna, std::unique_ptr<structure::Structure>* s = nullptr) {
  SetRna(rna);
  return Fold(s);
}

energy_t FoldBruteForce(std::unique_ptr<structure::Structure>* s = nullptr);

inline energy_t FoldBruteForce(const rna_t& rna, std::unique_ptr<structure::Structure>* s = nullptr) {
  SetRna(rna);
  return FoldBruteForce(s);
}

}
}

#endif //MEMERNA_FOLD_H
