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
// Access the array g_terminal[G][X][Y][U] since G, X, Y, U would occur in that order in memory.
extern energy_t g_stack[4][4][4][4];

// Terminal mismatch:
extern energy_t g_terminal[4][4][4][4];

// Internal loop related:
extern energy_t g_internal_init[INITIATION_CACHE_SZ];
extern energy_t g_internal_1x1[4][4][4][4][4][4];
extern energy_t g_internal_1x2[4][4][4][4][4][4][4];
extern energy_t g_internal_2x2[4][4][4][4][4][4][4][4];
extern energy_t g_internal_2x3_mismatch[4][4][4][4];
extern energy_t g_internal_other_mismatch[4][4][4][4];
extern energy_t g_internal_asym;
extern energy_t g_internal_augu_penalty;

// Bulge loop related:
extern energy_t g_bulge_init[INITIATION_CACHE_SZ];
extern energy_t g_bulge_special_c;

// Hairpin loop related:
extern energy_t g_hairpin_init[INITIATION_CACHE_SZ];
extern energy_t g_hairpin_uu_ga_first_mismatch, g_hairpin_gg_first_mismatch,
    g_hairpin_special_gu_closure, g_hairpin_c3_loop, g_hairpin_all_c_a, g_hairpin_all_c_b;
extern std::unordered_map<std::string, energy_t> g_hairpin;

// Multiloop hack model:
extern energy_t g_multiloop_hack_a, g_multiloop_hack_b;

// Dangles:
// X, G, U
extern energy_t g_dangle5[4][4][4];
// G, U, X
extern energy_t g_dangle3[4][4][4];

// Coaxial stacking:
extern energy_t g_coax_mismatch_non_contiguous, g_coax_mismatch_wc_bonus, g_coax_mismatch_gu_bonus;

// AU/GU penalty
extern energy_t g_augu_penalty;

// Global data variables.
extern primary_t r;

inline void SetPrimary(const primary_t& primary) { r = primary; }

}

#endif  //MEMERNA_GLOBALS_H
