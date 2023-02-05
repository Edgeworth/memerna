// Copyright 2016 Eliot Courtney.
#ifndef COMMON_TEST_H_
#define COMMON_TEST_H_

#include <cstdint>
#include <string>

#include "api/energy/model.h"
#include "models/t04/energy/model.h"
#include "models/t22/energy/model.h"

namespace mrna {

inline constexpr uint32_t T04P1_MODEL_HASH = 0x443cf312;
inline constexpr uint32_t T04P2_MODEL_HASH = 0xdfdf6e87;
inline constexpr uint32_t T12P2_MODEL_HASH = 0xdbcc795b;
inline constexpr uint32_t T22P2_MODEL_HASH = 0xdbcc795b;

// Make sure to use Range(0, NUM_TEST_MODELS) if making a parameterised test
// with all models in test_ems, since test_ems is initialized at runtime.
inline constexpr int NUM_TEST_MODELS = 5;
inline erg::EnergyModelPtr test_ems[NUM_TEST_MODELS];
inline md::t04::Model::Ptr test_t04_ems[NUM_TEST_MODELS];
inline md::t22::Model::Ptr test_t22_ems[NUM_TEST_MODELS];

#if ENERGY_PRECISION == 1
inline md::t04::Model::Ptr t04p1;
#elif ENERGY_PRECISION == 2
inline md::t04::Model::Ptr t04p2;
inline md::t04::Model::Ptr t12p2;
inline md::t22::Model::Ptr t22p2;

// NEWMODEL: Add here.

#endif

void InitTest(const std::string& data_dir);

}  // namespace mrna

#endif  // COMMON_TEST_H_
