// Copyright 2022 E.
#include "common_test.h"

namespace mrna {

void InitTest(const std::string& data_dir) {
  int last_test_ems_idx = 0;
  int last_test_t04_ems_idx = 0;

  // NEWMODEL: Add here.
#if ENERGY_PRECISION == 1
  mrna::t04p1 = mrna::erg::t04::Model::FromDir(mrna::erg::ModelPath(data_dir, "t04p1"));
  verify(mrna::t04p1->Checksum() == mrna::T04P1_MODEL_HASH, "Expected t04p1 model hash {}, got {}",
      mrna::T04P1_MODEL_HASH, mrna::t04p1->Checksum());

  mrna::test_ems[last_test_ems_idx++] = mrna::t04p1;
  mrna::test_t04_ems[last_test_t04_ems_idx++] = mrna::t04p1;

#elif ENERGY_PRECISION == 2
  mrna::t04p2 = mrna::erg::t04::Model::FromDir(mrna::erg::ModelPath(data_dir, "t04p2"));
  verify(mrna::t04p2->Checksum() == mrna::T04P2_MODEL_HASH, "Expected t04p2 model hash {}, got {}",
      mrna::T04P2_MODEL_HASH, mrna::t04p2->Checksum());

  mrna::t12p2 = mrna::erg::t04::Model::FromDir(mrna::erg::ModelPath(data_dir, "t12p2"));
  verify(mrna::t12p2->Checksum() == mrna::T12P2_MODEL_HASH, "Expected t12p2 model hash {}, got {}",
      mrna::T12P2_MODEL_HASH, mrna::t12p2->Checksum());

  // TODO(0): uncomment.
  // mrna::t22p2 = mrna::erg::t22::Model::FromDir(mrna::erg::ModelPath(data_dir, "t22p2"));
  // verify(mrna::t22p2->Checksum() == mrna::T22P2_MODEL_HASH, "Expected t22 model hash {}, got {}",
  //     mrna::T22P2_MODEL_HASH, mrna::t22p2->Checksum());

  mrna::test_ems[last_test_ems_idx++] = mrna::t04p2;
  mrna::test_ems[last_test_ems_idx++] = mrna::t12p2;
  mrna::test_ems[last_test_ems_idx++] = mrna::t22p2;
  mrna::test_t04_ems[last_test_t04_ems_idx++] = mrna::t04p2;

#endif

  for (int i = last_test_ems_idx; i < mrna::NUM_TEST_MODELS; ++i)
    mrna::test_ems[i] = mrna::erg::Random(mrna::erg::ModelKind::T04_LIKE, i);

  for (int i = last_test_t04_ems_idx; i < mrna::NUM_TEST_MODELS; ++i)
    mrna::test_t04_ems[i] = mrna::erg::t04::Model::Random(i);
}

}  // namespace mrna
