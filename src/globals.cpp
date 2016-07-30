#include "globals.h"

namespace memerna {

energy_t stacking_e[4][4][4][4];

energy_t terminal_e[4][4][4][4];

energy_t internal_init[INITIATION_CACHE_SZ];
energy_t internal_1x1[4][4][4][4][4][4];
energy_t internal_1x2[4][4][4][4][4][4][4];
energy_t internal_2x2[4][4][4][4][4][4][4][4];
energy_t internal_2x3_mismatch[4][4][4][4];
energy_t internal_other_mismatch[4][4][4][4];
energy_t internal_asym;
energy_t internal_augu_penalty;
energy_t internal_mismatch_1xk;

energy_t bulge_init[INITIATION_CACHE_SZ];
energy_t bulge_special_c;

energy_t hairpin_init[INITIATION_CACHE_SZ];
energy_t hairpin_uu_ga_first_mismatch, hairpin_gg_first_mismatch,
    hairpin_special_gu_closure, hairpin_c3_loop, hairpin_all_c_a, hairpin_all_c_b;
std::unordered_map<std::string, energy_t> hairpin_e;

energy_t multiloop_hack_a, multiloop_hack_b;

energy_t multiloop_t99_a, multiloop_t99_b, multiloop_t99_c;

energy_t dangle5_e[4][4][4];
energy_t dangle3_e[4][4][4];

energy_t coax_mismatch_non_contiguous, coax_mismatch_wc_bonus, coax_mismatch_gu_bonus;

rna_t r;
std::vector<int> p;

}
