#ifndef MEMERNA_GLOBALS_H
#define MEMERNA_GLOBALS_H

#include <string>
#include <vector>
#include <unordered_map>
#include "base.h"

namespace memerna {

// Globals.
// ----!!!!IF YOU UPDATE THESE UPDATE SerialiseEnergyModel() and LoadRandomEnergyModel() IN common.!!!!----
const int INITIATION_CACHE_SZ = 31;

// Stacking related:
// Note that the order of indices is always the order the RNA would be accessed in in memory.
// E.g. For a terminal mismatch:
// 5'-> G X 3'->
// <-3' U Y <-5'
// Access the array terminal_e[G][X][Y][U] since G, X, Y, U would occur in that order in memory.
extern energy_t stacking_e[4][4][4][4];

// Terminal mismatch:
extern energy_t terminal_e[4][4][4][4];

// Internal loop related:
extern energy_t internal_init[INITIATION_CACHE_SZ];
extern energy_t internal_1x1[4][4][4][4][4][4];
extern energy_t internal_1x2[4][4][4][4][4][4][4];
extern energy_t internal_2x2[4][4][4][4][4][4][4][4];
extern energy_t internal_2x3_mismatch[4][4][4][4];
extern energy_t internal_other_mismatch[4][4][4][4];
extern energy_t internal_asym;
extern energy_t internal_augu_penalty;
extern energy_t internal_mismatch_1xk;

// Bulge loop related:
extern energy_t bulge_init[INITIATION_CACHE_SZ];
extern energy_t bulge_special_c;

// Hairpin loop related:
extern energy_t hairpin_init[INITIATION_CACHE_SZ];
extern energy_t hairpin_uu_ga_first_mismatch, hairpin_gg_first_mismatch,
    hairpin_special_gu_closure, hairpin_c3_loop, hairpin_all_c_a, hairpin_all_c_b;
extern std::unordered_map<std::string, energy_t> hairpin_e;

// Multiloop hack model:
extern energy_t multiloop_hack_a, multiloop_hack_b;

// Multiloop T99 model:
extern energy_t multiloop_t99_a, multiloop_t99_b, multiloop_t99_c;

// Dangles:
// X, G, U
extern energy_t dangle5_e[4][4][4];
// G, U, X
extern energy_t dangle3_e[4][4][4];

// Coaxial stacking:
extern energy_t coax_mismatch_non_contiguous, coax_mismatch_wc_bonus, coax_mismatch_gu_bonus;

// AU/GU penalty
extern energy_t augu_penalty;

// Global data variables.
extern rna_t r;
extern std::vector<int> p;

inline void SetRna(const rna_t& rna) {
  r = rna;
}

inline void SetFoldedRna(const rna_t& rna, const std::vector<int>& pairs) {
  r = rna;
  p = pairs;
}

inline void SetFoldedRna(const folded_rna_t& frna) {
  SetFoldedRna(frna.r, frna.p);
}

}

#endif  //MEMERNA_GLOBALS_H
