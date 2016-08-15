#ifndef MEMERNA_FOLD_H
#define MEMERNA_FOLD_H

#include "common.h"
#include "energy/energy.h"

namespace memerna {
namespace fold {

energy_t Fold();
inline energy_t Fold(const rna_t& rna) {
  SetRna(rna);
  return Fold();
}
energy_t FoldBruteForce();

}
}

#endif //MEMERNA_FOLD_H
