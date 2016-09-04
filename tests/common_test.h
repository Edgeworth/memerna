#ifndef MEMERNA_COMMON_TEST_H
#define MEMERNA_COMMON_TEST_H

#include <cstdint>

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
const char* const ENERGY_MODEL_PATH = "data";

}

#endif //MEMERNA_COMMON_TEST_H
