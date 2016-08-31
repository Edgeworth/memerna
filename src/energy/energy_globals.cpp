#include "energy/energy_globals.h"

namespace memerna {
namespace energy {

energy_t g_stack[4][4][4][4];

energy_t g_terminal[4][4][4][4];

energy_t g_internal_init[INITIATION_CACHE_SZ];
energy_t g_internal_1x1[4][4][4][4][4][4];
energy_t g_internal_1x2[4][4][4][4][4][4][4];
energy_t g_internal_2x2[4][4][4][4][4][4][4][4];
energy_t g_internal_2x3_mismatch[4][4][4][4];
energy_t g_internal_other_mismatch[4][4][4][4];
energy_t g_internal_asym;
energy_t g_internal_augu_penalty;

energy_t g_bulge_init[INITIATION_CACHE_SZ];
energy_t g_bulge_special_c;

energy_t g_hairpin_init[INITIATION_CACHE_SZ];
energy_t g_hairpin_uu_ga_first_mismatch, g_hairpin_gg_first_mismatch,
    g_hairpin_special_gu_closure, g_hairpin_c3_loop, g_hairpin_all_c_a, g_hairpin_all_c_b;
std::unordered_map<std::string, energy_t> g_hairpin;

energy_t g_multiloop_hack_a, g_multiloop_hack_b;

energy_t g_dangle5[4][4][4];
energy_t g_dangle3[4][4][4];

energy_t g_coax_mismatch_non_contiguous, g_coax_mismatch_wc_bonus, g_coax_mismatch_gu_bonus;

energy_t g_augu_penalty;

}
}
