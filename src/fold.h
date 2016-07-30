#ifndef MEMERNA_FOLD_H
#define MEMERNA_FOLD_H

#include "common.h"
#include "energy.h"

namespace memerna {
namespace fold {

// This function is not re-entrant.
energy_t Fold(const rna_t& rna, std::unique_ptr<structure::Structure>* s = nullptr);

}
}

#endif //MEMERNA_FOLD_H
