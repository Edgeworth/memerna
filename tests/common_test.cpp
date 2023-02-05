// Copyright 2022 Eliot Courtney.
#include "common_test.h"

#include <variant>

#include "api/energy/energy_cfg.h"

namespace mrna {

void InitTest(const std::string& data_dir) {
  int last_test_ems_idx = 0;
  int last_test_t04_ems_idx = 0;
  int last_test_t22_ems_idx = 0;

  // NEWMODEL: Add here.
#if ENERGY_PRECISION == 1
  mrna::t04p1 = mrna::md::t04::Model::FromDir(mrna::erg::ModelPath(data_dir, "t04p1"));

  mrna::test_ems[last_test_ems_idx++] = mrna::t04p1;
  mrna::test_t04_ems[last_test_t04_ems_idx++] = mrna::t04p1;

#elif ENERGY_PRECISION == 2
  mrna::t04p2 = mrna::md::t04::Model::FromDir(mrna::erg::ModelPath(data_dir, "t04p2"));
  mrna::t12p2 = mrna::md::t04::Model::FromDir(mrna::erg::ModelPath(data_dir, "t12p2"));
  mrna::t22p2 = mrna::md::t22::Model::FromDir(mrna::erg::ModelPath(data_dir, "t22p2"));

  mrna::test_ems[last_test_ems_idx++] = mrna::t04p2;
  mrna::test_ems[last_test_ems_idx++] = mrna::t12p2;
  mrna::test_ems[last_test_ems_idx++] = mrna::t22p2;
  mrna::test_t04_ems[last_test_t04_ems_idx++] = mrna::t04p2;
  mrna::test_t04_ems[last_test_t04_ems_idx++] = mrna::t12p2;
  mrna::test_t22_ems[last_test_t22_ems_idx++] = mrna::t22p2;

#endif

  for (int i = last_test_ems_idx; i < mrna::NUM_TEST_MODELS; ++i)
    mrna::test_ems[i] = mrna::erg::Random(mrna::erg::ModelKind::T04_LIKE, i);

  for (int i = last_test_t04_ems_idx; i < mrna::NUM_TEST_MODELS; ++i)
    mrna::test_t04_ems[i] = mrna::md::t04::Model::Random(i);

  for (int i = last_test_t22_ems_idx; i < mrna::NUM_TEST_MODELS; ++i)
    mrna::test_t22_ems[i] = mrna::md::t22::Model::Random(i);
}

}  // namespace mrna
