#ifndef MEMERNA_GLOBALS_H
#define MEMERNA_GLOBALS_H

#include <string>
#include <vector>
#include <unordered_map>
#include "base.h"
#include "energy.h"

namespace memerna {

// Globals.
const int INITIATION_CACHE_SZ = 31;

// Stacking related:
// Note that the order of indices is always the order the RNA would be accessed in in memory.
// E.g. For a terminal mismatch:
// 5'-> G X 3'->
// <-3' U Y <-5'
// Access the array terminal_e[G][X][Y][U] since G, X, Y, U would occur in that order in memory.
extern energy::energy_t stacking_e[4][4][4][4];

// Terminal mismatch:
extern energy::energy_t terminal_e[4][4][4][4];

// Internal loop related:
extern energy::energy_t internal_init[INITIATION_CACHE_SZ];
extern energy::energy_t internal_1x1[4][4][4][4][4][4];
extern energy::energy_t internal_1x2[4][4][4][4][4][4][4];
extern energy::energy_t internal_2x2[4][4][4][4][4][4][4][4];
extern energy::energy_t internal_2x3_mismatch[4][4][4][4];
extern energy::energy_t internal_other_mismatch[4][4][4][4];
extern energy::energy_t internal_asym;
extern energy::energy_t internal_augu_penalty;
extern energy::energy_t internal_mismatch_1xk;

// Bulge loop related:
extern energy::energy_t bulge_init[INITIATION_CACHE_SZ];
extern energy::energy_t bulge_special_c;

// Hairpin loop related:
extern energy::energy_t hairpin_init[INITIATION_CACHE_SZ];
extern energy::energy_t hairpin_uu_ga_first_mismatch, hairpin_gg_first_mismatch,
    hairpin_special_gu_closure, hairpin_c3_loop, hairpin_all_c_a, hairpin_all_c_b;
extern std::unordered_map<std::string, energy::energy_t> hairpin_e;

// Multiloop related:
extern energy::energy_t multiloop_a, multiloop_b, multiloop_c;

// Dangles:
// X, G, U
extern energy::energy_t dangle5_e[4][4][4];
// G, U, X
extern energy::energy_t dangle3_e[4][4][4];

// Coaxial stacking:
extern energy::energy_t coax_mismatch_non_contiguous, coax_mismatch_wc_bonus, coax_mismatch_gu_bonus;

// Global data variables.
extern rna_t r;
extern std::vector<int> p;

}

#endif //MEMERNA_GLOBALS_H
