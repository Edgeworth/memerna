#include "globals.h"

namespace memerna {

energy::energy_t stacking_e[4][4][4][4];

energy::energy_t terminal_e[4][4][4][4];

energy::energy_t internal_init[INITIATION_CACHE_SZ];

energy::energy_t bulge_init[INITIATION_CACHE_SZ];

energy::energy_t hairpin_init[INITIATION_CACHE_SZ];
energy::energy_t hairpin_uu_ga_first_mismatch, hairpin_gg_first_mismatch,
    hairpin_special_gu_closure, hairpin_c3_loop, hairpin_all_c_a, hairpin_all_c_b;
std::unordered_map<std::string, energy::energy_t> hairpin_e;

rna_t r;
std::vector<int> p;

}
