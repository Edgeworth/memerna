// Copyright 2016 Eliot Courtney.
#ifndef TESTS_INIT_H_
#define TESTS_INIT_H_

#include <string>
#include <tuple>
#include <variant>
#include <vector>

#include "api/ctx/backend.h"
#include "api/ctx/backend_cfg.h"
#include "backends/stack/energy/model.h"
#include "model/primary.h"
#include "model/secondary.h"

namespace mrna {

// BackendCfg constants for IsSupported checks in tests.
// Uses monostate for data_src since we only need energy_model for checks.
inline const BackendCfg kT04Cfg{.energy_model = erg::EnergyModelKind::T04,
    .precision = ENERGY_PRECISION,
    .data_src = std::monostate{}};

inline const BackendCfg kT12Cfg{.energy_model = erg::EnergyModelKind::T12,
    .precision = ENERGY_PRECISION,
    .data_src = std::monostate{}};

inline const BackendCfg kT22Cfg{.energy_model = erg::EnergyModelKind::T22,
    .precision = ENERGY_PRECISION,
    .data_src = std::monostate{}};

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
// with all models in <backend>_ms, since they are initialized at runtime.
inline constexpr int NUM_TEST_MODELS = 5;

inline md::base::Model::Ptr base_t04;
inline md::base::opt::Model::Ptr baseopt_t04;
inline md::stack::Model::Ptr stack_t04;

inline std::vector<md::base::Model::Ptr> base_ms;
inline std::vector<md::base::opt::Model::Ptr> baseopt_ms;
inline std::vector<md::stack::Model::Ptr> stack_ms;
// NEWBACKEND: Add here.

inline constexpr int NUM_T04_MODELS = 3;
inline std::vector<BackendModelPtr> t04_ms;

inline constexpr int NUM_T12_MODELS = 3;
inline std::vector<BackendModelPtr> t12_ms;

inline constexpr int NUM_T22_MODELS = 1;
inline std::vector<BackendModelPtr> t22_ms;
// NEWMODEL: Add here.

void InitTest(const std::string& data_dir);

}  // namespace mrna

#endif  // TESTS_INIT_H_
