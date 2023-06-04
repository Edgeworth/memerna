// Copyright 2016 Eliot Courtney.
#ifndef COMMON_TEST_H_
#define COMMON_TEST_H_

#include <string>
#include <tuple>

#include "api/energy/model.h"
#include "model/primary.h"
#include "model/secondary.h"
#include "models/t04/energy/model.h"
#include "models/t22/energy/model.h"

namespace mrna {

#define EXPECT_REL_EQ(a, b)                                                    \
  do {                                                                         \
    auto acopy = (a);                                                          \
    auto bcopy = (b);                                                          \
    EXPECT_TRUE(rel_eq(acopy, bcopy))                                          \
        << std::setprecision(FLOAT_PRECISION + 1) << acopy << " != " << bcopy; \
  } while (0)

extern std::tuple<Primary, Secondary> kNNDBHairpin1;
extern std::tuple<Primary, Secondary> kNNDBHairpin2;
extern std::tuple<Primary, Secondary> kNNDBHairpin3;
extern std::tuple<Primary, Secondary> kNNDBHairpin4;
extern std::tuple<Primary, Secondary> kNNDBHairpin5;
extern std::tuple<Primary, Secondary> kNNDBBulge1;
extern std::tuple<Primary, Secondary> kNNDBBulge2;
extern std::tuple<Primary, Secondary> kNNDBInternal2x3;
extern std::tuple<Primary, Secondary> kNNDBInternal1x5;
extern std::tuple<Primary, Secondary> kNNDBInternal2x2;
extern std::tuple<Primary, Secondary> kBulge1;
extern std::tuple<Primary, Secondary> kInternal1;
extern std::tuple<Primary, Secondary> k16sHSapiens3;

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
