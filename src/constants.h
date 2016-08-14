#ifndef MEMERNA_CONSTANTS_H
#define MEMERNA_CONSTANTS_H

#include "common.h"

namespace memerna {
namespace constants {

// Don't change this value. Plays nice with memset.
const energy_t MAX_E = 0x0F0F0F0F;
const energy_t CAP_E = 0x07070707;

// -----------------------------------------------
// Values affecting the energy model:
const int HAIRPIN_MIN_SZ = 3;
// N.B. This is for kcal/mol so it's not 8.315.
const double R = 1.985877534e-3;
// This is 37 degrees Celsius. Changing this is not a good idea.
const double T = 310.15;
// Ninio maximum asymmetry.
const energy_t NINIO_MAX_ASYM = 30;
// Maximum size of a twoloop.
const int TWOLOOP_MAX_SZ = 30;

}
}

#endif //MEMERNA_CONSTANTS_H
