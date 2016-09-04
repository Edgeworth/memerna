#ifndef MEMERNA_COMMON_TEST_H
#define MEMERNA_COMMON_TEST_H

#include <cstdint>
#include "energy/energy_model.h"

#define ONLY_FOR_THIS_MODEL(em_, hash_) \
  do { \
    auto our_hash = em_.Checksum(); \
    if (our_hash != hash_) { \
      printf("Skipping energy model specific tests: %#010x != " #hash_ " (%#010x).\n", our_hash, hash_); \
      return; \
    } \
  } while(0)

namespace memerna {

const uint32_t T04_MODEL_HASH = 0x03b94db8;
extern energy::EnergyModel g_em;
extern std::vector<energy::EnergyModel> g_ems;

}

#endif //MEMERNA_COMMON_TEST_H
