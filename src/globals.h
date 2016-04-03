#ifndef MEMERNA_GLOBALS_H
#define MEMERNA_GLOBALS_H

#include <string>
#include <vector>
#include <unordered_map>
#include "base.h"
#include "energy.h"

namespace memerna {

// Globals.
// A X Y A
const int INITIATION_CACHE_SZ = 31;

// Stacking related:
extern energy::energy_t stacking_e[4][4][4][4];

extern energy::energy_t terminal_e[4][4][4][4];

// Internal loop related:
extern energy::energy_t internal_init[INITIATION_CACHE_SZ];
extern energy::energy_t internal_1x1[4][4][4][4][4][4];
extern energy::energy_t internal_1x2[4][4][4][4][4][4][4];
extern energy::energy_t internal_2x2[4][4][4][4][4][4][4][4];

// Bulge loop related:
extern energy::energy_t bulge_init[INITIATION_CACHE_SZ];
extern energy::energy_t bulge_special_c;

// Hairpin loop related:
extern energy::energy_t hairpin_init[INITIATION_CACHE_SZ];
extern energy::energy_t hairpin_uu_ga_first_mismatch, hairpin_gg_first_mismatch,
    hairpin_special_gu_closure, hairpin_c3_loop, hairpin_all_c_a, hairpin_all_c_b;
extern std::unordered_map<std::string, energy::energy_t> hairpin_e;

// Global data variables.
extern rna_t r;
extern std::vector<int> p;

}

#endif //MEMERNA_GLOBALS_H
