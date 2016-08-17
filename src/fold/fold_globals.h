#ifndef MEMERNA_FOLD_GLOBALS_H
#define MEMERNA_FOLD_GLOBALS_H

#include "common.h"

namespace memerna {
namespace fold {

// For checks to make sure everything has been initialised.
extern bool g_fold_init;
extern energy_t g_augubranch[4][4];

}
}

#endif //MEMERNA_FOLD_GLOBALS_H
