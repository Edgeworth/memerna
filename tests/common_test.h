#ifndef MEMERNA_COMMON_TEST_H
#define MEMERNA_COMMON_TEST_H

#include <cstdint>

#define ONLY_FOR_THIS_MODEL(hash) \
  do { \
    auto our_hash = energy::EnergyModelChecksum(); \
    if (our_hash != hash) { \
      printf("Skipping energy model specific tests: %#010x != " #hash " (%#010x).\n", our_hash, hash); \
      return; \
    } \
  } while(0)

namespace memerna {

const uint32_t T04_MODEL_HASH = 0x03b94db8;

}

#endif //MEMERNA_COMMON_TEST_H
