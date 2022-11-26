// Copyright 2022 Eliot Courtney.
#include "common_test.h"

namespace mrna {

void InitTest(const std::string& data_dir) {
#if ENERGY_PRECISION == 1

  mrna::t04_p1 = mrna::erg::t04::Model::FromDir(mrna::erg::ModelPath(data_dir, "t04_p1"));
  verify(mrna::t04_p1->Checksum() == mrna::T04_P1_MODEL_HASH,
      "Expected t04_p1 model hash %d, got %d", mrna::T04_P1_MODEL_HASH, mrna::t04_p1->Checksum());

  mrna::test_ems[0] = mrna::t04_p1;
  mrna::test_t04_ems[0] = mrna::t04_p1;

#elif ENERGY_PRECISION == 2

  // TODO(0): support t12
  mrna::t04_p2 = mrna::erg::t04::Model::FromDir(mrna::erg::ModelPath(data_dir, "t04_p2"));
  verify(mrna::t04_p2->Checksum() == mrna::T04_P2_MODEL_HASH,
      "Expected t04_p2 model hash %d, got %d", mrna::T04_P2_MODEL_HASH, mrna::t04_p2->Checksum());

  mrna::test_ems[0] = mrna::t04_p2;
  mrna::test_t04_ems[0] = mrna::t04_p2;

#endif

  // TODO(1): support t22
  for (int i = 1; i < mrna::NUM_TEST_MODELS; ++i) {
    mrna::test_ems[i] = mrna::erg::Random(mrna::erg::ModelKind::T04_LIKE, i);
    mrna::test_t04_ems[i] = mrna::erg::t04::Model::Random(i);
  }
}

}  // namespace mrna
